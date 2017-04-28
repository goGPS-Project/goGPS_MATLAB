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
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.1 beta
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
    
    % Default values for each field - useful to restore corrupted field
    properties (Constant, Access = 'protected')
        
        % PROJECT                
        PRJ_NAME = 'Defauld DD';  % Name of the project
        PRJ_HOME = [IO_Settings.DEFAULT_DIR_IN 'project' filesep 'default_DD' filesep]; % Location of the project <relative path from goGPS folder> 
        CUR_INI = [IO_Settings.DEFAULT_DIR_IN 'project' filesep 'default_DD' filesep 'Config' filesep 'settings.ini']; % Location of the current ini file
        
        % DEPRECATE
        INPUT_FILE_INI_PATH = ''; % deprecate INI - it contains some additional setting (yet not imported in the new settings system)

        % RECEIVERS        
        SSS_DATE_START = GPS_Time(); % Start of the processing session
        SSS_DATE_STOP = GPS_Time();  % End of the processing session
        SSS_ID_LIST = '0';   % id character sequence to be use for the session $(S) special keyword
        SSS_ID_START = '0';  % first session id (char of sss_id_list)
        SSS_ID_STOP = '0';   % last session id (char of sss_id_list)
        
        OBS_DIR = '../data/project/default_DD/RINEX/' % 
        OBS_NAME = {'yamatogawa_master.obs' 'yamatogawa_rover.obs'} ;
        REC_TARGET = 0;
        REC_MASTER = 1;
        REC_REFERENCE = 2;
        OBS_TYPE_LIST = {'Target', ...
                         'Master', ...
                         'SEID reference'};
        OBS_TYPE = [1 0];

        ATX_DIR = [IO_Settings.DEFAULT_DIR_IN 'antenna' filesep 'ATX' filesep]; % Location of the antex files
        ATX_NAME = 'I08.ATX';    % Name antex file
        
        XYZ_ANT = zeros(3, 1);
        XYZ_EV_POINT = zeros(3, 1);

        % COMPUTATION CENTERS
        % With official products for orbits and clocks
        PREFERRED_ARCHIVE = {'cddis', 'igscb', 'custom'}
        PREFERRED_GPS = {'igs', 'emx', 'gfz'}
        PREFERRED_GLO = {'igs', 'emx', 'gfz'}
        PREFERRED_MXD = {'gbm'}
        PREFERRED_CLK = {'clk_30s', 'clk'}
        PREFERRED_EPH = {'final', 'rapid', 'ultra', 'broadcast'}
        
        CUSTOM_ADDR = 'cddis.gsfc.nasa.gov/'
        CUSTOM_PORT = '21'
        CUSTOM_PATH = 'pub/gps/products/'
        CUSTOM_NAME_EPH = '${WWWW}/igs${WWWWD}.sp3';
        CUSTOM_NAME_ERP = '${WWWW}/igs${WWWWD}.erp';
        CUSTOM_NAME_CLK = '${WWWW}/igs${WWWWD}.clk_30s';                
                
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
        GEOID_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'geoid' filesep]; % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs        
        
        % DTM (SET PATH AND LOAD PARAMETER FILES)        
        DTM_DIR = [IO_Settings.DEFAULT_DIR_IN 'reference' filesep 'DTM' filesep]; % Path to DTM folder containing DTM files
        % UI IMAGES
        IMG_DIR = [IO_Settings.DEFAULT_DIR_IN 'img' filesep];  % Path to images used by the interface
        IMG_LOGO64 = [IO_Settings.IMG_DIR 'goGPS_logo_64.png'];     % Location of the goGPS logo 64x64
        OUT_DIR = [IO_Settings.DEFAULT_DIR_OUT  'project' filesep 'default_DD' filesep 'out' filesep]; % Directory containing the output of the project
        OUT_PREFIX = 'out';  % Every time a solution is computed a folder with prefix followed by the run number is created
        RUN_COUNTER = 0;     % This parameter store the current run number
        
        % EXTERNAL INFO as imported from the input ini file does not have default values
    end
    
    properties (Constant, Access = 'protected')
        % id to string of out modes
        DEFAULT_DIR_IN = ['..' filesep 'data' filesep];
        DEFAULT_DIR_OUT = ['..' filesep 'data' filesep];  
    end
    
    properties (Constant, Access = 'public')
        % Location of the latest project (the ini contains just a reference to the default project ini file - that is actually a settings file
        LAST_SETTINGS = [IO_Settings.DEFAULT_DIR_IN 'last_settings.ini'];
    end
    

    properties (SetAccess = protected, GetAccess = protected)
        % Location of the current ini file
        cur_ini = IO_Settings.CUR_INI;
    end
    
    properties (SetAccess = protected, GetAccess = public)
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
        preferred_clk = IO_Settings.PREFERRED_CLK;          % kind of orbits to prefer
        
        % Custom entry for a server
        
        custom_addr = IO_Settings.CUSTOM_ADDR; % ftp address for a custom server
        custom_port = IO_Settings.CUSTOM_PORT; % ftp port for a custom server
        custom_path = IO_Settings.CUSTOM_PATH; % remote path 
        custom_name_eph = IO_Settings.CUSTOM_NAME_EPH; % name of the ephemeris file in the remote dir
        custom_name_erp = IO_Settings.CUSTOM_NAME_ERP; % name of the Earth Rotation Parameters file in the remote dir
        custom_name_clk = IO_Settings.CUSTOM_NAME_CLK; % name of the clock file in the remote dir

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
        
        % Path to stations metereological files
        met_dir = IO_Settings.MET_DIR;
        % Location of the meteorological file
        met_name =  IO_Settings.MET_NAME;

        % Path to stations ocean loading files
        ocean_dir = IO_Settings.OCEAN_DIR;

        %------------------------------------------------------------------
        % REFERENCE
        %------------------------------------------------------------------
        
        % Path to file containing the reference path
        ref_graph_file = IO_Settings.REF_GRAPH_FILE;

        % Path to ERP folder containing Earth Rotation Parameters (tipically realesed together with the orbits)
        erp_dir = IO_Settings.GEOID_DIR;
        
        % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        geoid_dir = IO_Settings.GEOID_DIR;
        
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
        
        % Location of the goGPS logo 64x64
        img_logo64 = IO_Settings.IMG_LOGO64;      
        
        %------------------------------------------------------------------
        % OUTPUT
        %------------------------------------------------------------------
        
        % Directory containing the output of the project
        out_dir = IO_Settings.OUT_DIR;
        
        % Every time a solution is computed a folder with prefix followed by the run number is created
        out_prefix = IO_Settings.OUT_PREFIX;
        
        % This parameter store the current run number
        run_counter = IO_Settings.RUN_COUNTER;
        
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
            this.img_logo64 = [this.img_dir 'goGPS_logo_64.png'];
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
                this.prj_home   = fnp.checkPath(settings.getData('prj_home'));                
                % COMPUTATION CENTERS
                this.preferred_archive = fnp.checkPath(settings.getData('preferred_archive'));
                this.preferred_gps = fnp.checkPath(settings.getData('preferred_gps'));
                this.preferred_glo = fnp.checkPath(settings.getData('preferred_glo'));
                this.preferred_mxd = fnp.checkPath(settings.getData('preferred_mxd'));
                this.preferred_eph = fnp.checkPath(settings.getData('preferred_eph'));                
                this.preferred_clk = fnp.checkPath(settings.getData('preferred_clk'));
                % Custom entry for a server                
                this.custom_addr = fnp.checkPath(settings.getData('custom_addr'));
                this.custom_port = fnp.checkPath(settings.getData('custom_port'));
                this.custom_path = fnp.checkPath(settings.getData('custom_path'));
                this.custom_name_eph = fnp.checkPath(settings.getData('custom_name_eph'));
                this.custom_name_erp = fnp.checkPath(settings.getData('custom_name_erp'));
                this.custom_name_clk = fnp.checkPath(settings.getData('custom_name_clk'));                
                % DEPRECATE
                this.input_file_ini_path = fnp.checkPath(settings.getData('input_file_ini_path'));
                % RECEIVERS
                this.sss_date_start = GPS_Time(datenum(settings.getData('sss_date_start')));
                this.sss_date_stop = GPS_Time(datenum(settings.getData('sss_date_stop')));
                this.sss_id_list = settings.getData('sss_id_list');
                this.sss_id_start = settings.getData('sss_id_start');
                this.sss_id_stop = settings.getData('sss_id_stop');                
                this.obs_dir = fnp.checkPath(settings.getData('obs_dir'));
                this.obs_name = fnp.checkPath(settings.getData('obs_name'));
                this.obs_type = settings.getData('obs_type');
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
                this.eph_dir    = fnp.checkPath(settings.getData('eph_dir'));
                this.eph_name   = fnp.checkPath(settings.getData('eph_name'));
                this.clk_dir    = fnp.checkPath(settings.getData('clk_dir'));
                this.clk_name   = fnp.checkPath(settings.getData('clk_name'));
                this.crx_dir    = fnp.checkPath(settings.getData('crx_dir'));
                this.dcb_dir    = fnp.checkPath(settings.getData('dcb_dir'));
                this.ems_dir    = fnp.checkPath(settings.getData('ems_dir'));
                % STATIONS
                this.crd_dir    = fnp.checkPath(settings.getData('crd_dir'));
                this.crd_name    = fnp.checkPath(settings.getData('crd_name'));
                this.met_dir    = fnp.checkPath(settings.getData('met_dir'));
                this.met_name    = fnp.checkPath(settings.getData('met_name'));
                this.ocean_dir  = fnp.checkPath(settings.getData('ocean_dir'));
                this.ocean_name  = fnp.checkPath(settings.getData('ocean_name'));
                % REFERENCE
                this.ref_graph_file  = fnp.checkPath(settings.getData('ref_graph_file'));
                this.geoid_dir  = fnp.checkPath(settings.getData('geoid_dir'));
                this.dtm_dir    = fnp.checkPath(settings.getData('dtm_dir'));
                % UI
                this.img_dir    = fnp.checkPath(settings.getData('img_dir'));
                this.img_logo64 = fnp.checkPath(settings.getData('img_logo64'));
                % OUTPUT
                this.out_dir = fnp.checkPath(settings.getData('out_dir'));
                this.out_prefix = fnp.checkPath(settings.getData('out_prefix'));
                this.run_counter = settings.getData('run_counter');
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
                this.preferred_clk = settings.preferred_clk;
                % Custom entry for a server
                this.custom_addr = settings.custom_addr;
                this.custom_port = settings.custom_port;
                this.custom_path = settings.custom_path;
                this.custom_name_eph = settings.custom_name_eph;
                this.custom_name_erp = settings.custom_name_erp;
                this.custom_name_clk = settings.custom_name_clk;
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
                this.obs_type = settings.obs_type;
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
                this.ocean_dir   = settings.ocean_dir;
                this.ocean_name  = settings.ocean_name;
                % REFERENCE
                this.ref_graph_file = settings.ref_graph_file;
                this.geoid_dir  = settings.geoid_dir;
                this.dtm_dir    = settings.dtm_dir;
                % UI
                this.img_dir    = settings.img_dir;
                this.img_logo64 = settings.img_logo64;
                % OUTPUT
                this.out_dir = settings.out_dir;
                this.out_prefix = settings.out_prefix;
                this.run_counter = settings.run_counter;
            end
            this.check();  
            this.updateExternals();
            this.eph_full_name = '';
            this.clk_full_name = '';
            this.updateObsFileName();
        end
        
        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- PROJECT --------------------------------------------------------------' 10 10];
            str = [str sprintf(' Project name:                                     %s\n', this.prj_name)];
            str = [str sprintf(' Project home:                                     %s\n', this.prj_home)];
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
            str = [str sprintf(' Directory of the observation files                %s \n', this.obs_dir)];
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
            str = [str sprintf(' Preferred order for clock products:               %s\n\n', strCell2Str(this.preferred_clk))];            
            str = [str sprintf(' Custom server parameters:\n')];
            str = [str sprintf('  address:                                         %s\n', this.custom_addr)];            
            str = [str sprintf('  port:                                            %s\n', this.custom_port)];            
            str = [str sprintf('  path:                                            %s\n', this.custom_path)];    
            str = [str sprintf('  eph name:                                        %s\n', this.custom_name_eph)];
            str = [str sprintf('  erp name:                                        %s\n', this.custom_name_erp)];
            str = [str sprintf('  clk name:                                        %s\n\n', this.custom_name_clk)];
            str = [str '---- INPUT FOLDERS: SATELLITE ---------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Ephemeris files:                     %s\n', this.eph_dir)];
            str = [str sprintf(' Name of Ephemeris files:                          %s\n', this.eph_name)];
            str = [str sprintf(' Directory of Satellite clock offsets:             %s\n', this.clk_dir)];
            str = [str sprintf(' Name of Satellite clock offsets:                  %s\n', this.clk_name)];
            str = [str sprintf(' Directory of CRX (satellite problems):            %s\n', this.crx_dir)];
            str = [str sprintf(' Directory of DCB (Differential Code Biases):      %s\n', this.dcb_dir)];
            str = [str sprintf(' Directory of EMS (EGNOS Message Server):          %s\n\n', this.ems_dir)];
            str = [str '---- INPUT FOLDERS: STATIONS ----------------------------------------------' 10 10];
            str = [str sprintf(' Directory of coordinates file:                    %s\n', this.crd_dir)];
            str = [str sprintf(' Name of coordinate (CRD) file:                    %s\n', this.crd_name)];
            str = [str sprintf(' Directory of meteorological data:                 %s\n', this.met_dir)];
            str = [str sprintf(' Name of meteorological (met) files:               %s\n', strCell2Str(this.met_name))];
            str = [str sprintf(' Directory of ocean loading files:                 %s\n', this.ocean_dir)];
            str = [str sprintf(' Name of ocean loading file:                       %s\n\n', this.ocean_name)];
            str = [str '---- INPUT FOLDERS: REFERENCE ---------------------------------------------' 10 10];
            str = [str sprintf(' File contraining the reference graph:             %s\n', this.ref_graph_file)];
            str = [str sprintf(' Directory of Geoid models:                        %s\n', this.geoid_dir)];
            str = [str sprintf(' Directory of DTM data:                            %s\n\n', this.dtm_dir)];
            str = [str '---- INPUT FOLDERS: UI ----------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of images for UI:                       %s\n', this.img_dir)];
            str = [str sprintf('  - Image of the goGPS logo:                       %s\n\n', this.img_logo64)];
            str = [str '---- OUTPUT SETTINGS ------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory containing the output of the project:   %s\n', this.out_dir)];
            str = [str sprintf(' Prefix of each run:                               %s\n', this.out_prefix)];
            str = [str sprintf(' Run counter:                                      %d\n\n', this.run_counter)];
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
            % PROJECT
            str_cell = Ini_Manager.toIniStringSection('PROJECT', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of the project', str_cell);
            str_cell = Ini_Manager.toIniString('prj_name', this.prj_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Home of the project', str_cell);
            str_cell = Ini_Manager.toIniString('prj_home', this.prj_home, str_cell);
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
            str_cell = Ini_Manager.toIniString('obs_dir', this.obs_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('File name of the receivers (can contain special keywords)', str_cell);
            str_cell = Ini_Manager.toIniString('obs_name', this.obs_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment(sprintf('Array of observations type (%s)', strCell2EnumStr(this.OBS_TYPE_LIST,', ')), str_cell);
            str_cell = Ini_Manager.toIniString('obs_type', this.obs_type, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of PCO - PCV antex (ATX) files', str_cell);
            str_cell = Ini_Manager.toIniString('atx_dir', this.atx_dir, str_cell);
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
            str_cell = Ini_Manager.toIniString('custom_name_erp', this.custom_name_erp, str_cell);
            str_cell = Ini_Manager.toIniString('custom_name_clk', this.custom_name_clk, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % SATELLITES
            str_cell = Ini_Manager.toIniStringSection('INPUT_SATELLITE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Ephemeris files', str_cell);
            str_cell = Ini_Manager.toIniString('eph_dir', this.eph_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of Ephemeris files - special keywords can be used', str_cell);
            str_cell = Ini_Manager.toIniStringComment('If not found, goGPS will try to download them following COMPUTATION_CENTER section', str_cell);
            str_cell = Ini_Manager.toIniString('eph_name', this.eph_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of clock offset files', str_cell);
            str_cell = Ini_Manager.toIniString('clk_dir', this.clk_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('If not found, goGPS will try to download them following COMPUTATION_CENTER section', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of clock offset files - special keywords can be used', str_cell);
            str_cell = Ini_Manager.toIniString('clk_name', this.clk_name, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of CRX files (containing satellite problems)', str_cell);
            str_cell = Ini_Manager.toIniString('crx_dir', this.crx_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DCB files (Differential Code Biases)', str_cell);
            str_cell = Ini_Manager.toIniString('dcb_dir', this.dcb_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of EMS files (EGNOS Message Server).', str_cell);
            str_cell = Ini_Manager.toIniString('ems_dir', this.ems_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % STATIONS
            str_cell = Ini_Manager.toIniStringSection('INPUT_STATIONS', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of coordinates files', str_cell);
            str_cell = Ini_Manager.toIniString('crd_dir', this.crd_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of coordinates (CRD) file', str_cell);
            str_cell = Ini_Manager.toIniString('crd_name', this.crd_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of metereological data', str_cell);
            str_cell = Ini_Manager.toIniString('met_dir', this.met_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Meteorological file', str_cell);
            str_cell = Ini_Manager.toIniString('met_name', this.met_name, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of ocean loading files', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_dir', this.ocean_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Name of ocean loading file', str_cell);
            str_cell = Ini_Manager.toIniString('ocean_name', this.ocean_name, str_cell);           
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % REFERENCE
            str_cell = Ini_Manager.toIniStringSection('INPUT_REFERENCE', str_cell);
            str_cell = Ini_Manager.toIniStringComment('File containing a graph of path constraints', str_cell);
            str_cell = Ini_Manager.toIniString('ref_graph_file', this.ref_graph_file, str_cell);            
            str_cell = Ini_Manager.toIniStringComment('Directory of Geoid files', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_dir', this.geoid_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DTM data', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', this.dtm_dir, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % UI
            str_cell = Ini_Manager.toIniStringSection('INPUT_UI', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of images for UI', str_cell);
            str_cell = Ini_Manager.toIniString('img_dir', this.img_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Path to the image of the logo 64x64 px', str_cell);
            str_cell = Ini_Manager.toIniString('img_logo64', this.img_logo64, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            % OUTPUT
            str_cell = Ini_Manager.toIniStringSection('OUTPUT', str_cell);

            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
            str_cell = Ini_Manager.toIniString('out_dir', this.out_dir, str_cell);
            str_cell = Ini_Manager.toIniString('out_prefix', this.out_prefix, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Current run number', str_cell);
            str_cell = Ini_Manager.toIniString('run_counter', this.run_counter, str_cell);
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);   
            
%             % EXTERNAL INFO
%             % this operation should not be blocking -> use try catch
%             % update external ini
%             try
%                 this.ext_ini = Ini_Manager(this.input_file_ini_path);
%                 
%                 str_cell = Ini_Manager.toIniStringSection('EXTERNAL FILE INPUT INI', str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('The information of this section are here present as imported from the external input ini file', str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('All the parameters here listed cannot be modified for import', str_cell);
%                 str_cell = Ini_Manager.toIniStringNewLine(str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('Variometric approach parameter', str_cell);
%                 str_cell = Ini_Manager.toIniString('variometric_time_step', this.ext_ini.getData('Variometric', 'timeStep'), str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('Navigational files', str_cell);
%                 str_cell = Ini_Manager.toIniString('nav_path', this.ext_ini.getData('Navigational', 'data_path'), str_cell);
%                 str_cell = Ini_Manager.toIniString('nav_file', this.ext_ini.getData('Navigational', 'file_name'), str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('Master/Target(if SEID) file', str_cell);
%                 str_cell = Ini_Manager.toIniString('master_target_path', this.ext_ini.getData('Master', 'data_path'), str_cell);
%                 str_cell = Ini_Manager.toIniString('master_target_file', this.ext_ini.getData('Master', 'file_name'), str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('Receivers/Source(if SEID) files', str_cell);
%                 str_cell = Ini_Manager.toIniString('receiver_source_number', this.ext_ini.getData('Receivers', 'nRec'), str_cell);
%                 str_cell = Ini_Manager.toIniString('receiver_source_path', this.ext_ini.getData('Receivers', 'data_path'), str_cell);
%                 str_cell = Ini_Manager.toIniString('receiver_source_file', this.ext_ini.getData('Receivers', 'file_name'), str_cell);
%                 str_cell = Ini_Manager.toIniStringComment('Binary files', str_cell);
%                 str_cell = Ini_Manager.toIniString('bin_path', this.ext_ini.getData('Bin', 'data_path'), str_cell);
%                 str_cell = Ini_Manager.toIniString('bin_file', this.ext_ini.getData('Bin', 'file_name'), str_cell);
%             catch ex
%                 this.logger.addWarning(sprintf('Exporting of external INFO failed - %s', ex.message));
%             end
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
                this.logger.addWarning(['Legacy import "IO file / folders" failed - ', ex.message])
            end
            this.check();
            this.updateExternals();
        end
    end
    
    % =========================================================================
    %  GETTERS
    % =========================================================================
    methods
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
        
        function dir = getLogo(this)
            % Get the logo of goGPS to be shown on the UI
            % SYNTAX: dir = this.getLogo()
            dir = this.img_logo64;
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
            % Get the file list of ephemeris files
            % SYNTAX: file_name = this.getFullNavClkPath(id)
            if isempty(this.clk_full_name)
                this.updateNavFileName();
            end
            file_name = this.clk_full_name;    
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
                        
        function out = getNavEphPath(this)
            % Get the path to the navigational files
            % SYNTAX: nav_path = this.getNavEphPath()
            out = File_Name_Processor.checkPath(fullfile(this.eph_dir, this.eph_name));
        end
        
        function out = getNavClkPath(this)
            % Get the path to the clock files
            % SYNTAX: nav_path = this.getNavClkPath()
            out = File_Name_Processor.checkPath(fullfile(this.clk_dir, this.clk_name));
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
        
        function out = getMetFile(this)
            % Get the path of the meteorological file
            % SYNTAX: file_path = this.getMetFile()
            if (isempty(this.met_name))
                out = '';
            else
                out = this.checkMetPath(strcat(this.met_dir, filesep, this.met_name));
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
        
        function out_prefix = getFullOutPrefix(this, date_start, date_stop, session_list, session_start, session_stop)
            % Get the full path of the outputs
            % SYNTAX: out_prefix = this.getFullOutPrefix()
            fnp = File_Name_Processor;
            out_prefix = fnp.checklPath(strcat(this.getOutDir, filesep, this.out_prefix));
            if nargin > 1
                out_prefix = fnp.dateKeyRepBatch(out_prefix, date_start, date_stop, session_list, session_start, session_stop);
            end            
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
                this.obs_full_name{i} = fnp.dateKeyRepBatch(fnp.checkPath(fullfile(this.obs_dir, this.obs_name{i})), this.sss_date_start,  this.sss_date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
            end
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

        function eph_full_name = getEphFileName(this, date_start, date_stop)
            % Get the full name of the ephemerides files (replacing special keywords)
            % SYNTAX: eph_full_name = getEphFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(fullfile(this.eph_dir, this.eph_name));
            step_sec = fnp.getStepSec(file_name);
            
            date_start = date_start.getCopy; date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
            date_stop = date_stop.getCopy; date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin            
            eph_full_name = fnp.dateKeyRepBatch(file_name, date_start,  date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
        
        function clk_full_name = getClkFileName(this, date_start, date_stop)
            % Get the full name of the clock offset files (replacing special keywords)
            % SYNTAX: clk_full_name = getClkFileName(this, date_start, date_stop)
            fnp = File_Name_Processor();
            file_name = fnp.checkPath(fullfile(this.clk_dir, this.clk_name));
            step_sec = fnp.getStepSec(file_name);
            
            date_start = date_start.getCopy; date_start.addIntSeconds(-step_sec); % Get navigational files with 6 hours of margin
            date_stop = date_stop.getCopy; date_stop.addIntSeconds(+step_sec); % Get navigational files with 6 hours of margin            
            clk_full_name = fnp.dateKeyRepBatch(file_name, date_start, date_stop, this.sss_id_list, this.sss_id_start, this.sss_id_stop);
        end
        
        
        function updateExternals(this)
            % Import the value of the external input files (stored in inputFile.ini)
            % SYNTAX: this.updateExternals();
            if ~isempty(this.input_file_ini_path)
                this.logger.addWarning('Legacy importing input ini files');
                this.ext_ini = Ini_Manager(this.input_file_ini_path);
                
                fnp = File_Name_Processor();
                
                % Import Rinex from old ini file
                dir_receiver = this.ext_ini.getData('Receivers', 'data_path');
                if ~isempty(dir_receiver)
                    this.obs_dir = dir_receiver;
                end
                dir_master = this.ext_ini.getData('Master', 'data_path');
                if ~(isempty(dir_master)) && ~strcmp(dir_master, this.obs_dir)
                    this.logger.addWarning('Importing legacy input file Master data_path seems different from Rover data_path - fix settings file manually');
                end
                
                % import Receivers/SEID names
                name_receiver = this.ext_ini.getData('Receivers', 'file_name');
                this.obs_name = {};
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
                this.updateObsFileName();
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
                this.prj_home = fnp.checkPath(path_parts{1:end-1}, filesep);
                this.prj_name = path_parts{end-1};
                this.logger.addMessage('Trying to guess project name / home / ini');
                this.logger.addMessage(sprintf(' name: %s', this.prj_name));
                this.logger.addMessage(sprintf(' home: %s', this.prj_home));
                this.logger.addMessage(sprintf(' ini:  %s', this.cur_ini));                
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

        function checkCellStringField(this, field_name, empty_is_valid, check_existence)
            % Check if a string field of the object is a valid string
            % To make the function works it is needed to have defined the default
            % value of the field as a constant with same name but upper case
            % SYNTAX: this.checkStringField(string_field_name, <empty_is_valid == false>, <check_existence == false>);            
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
            this.checkCellStringField('preferred_clk', false);
            
            this.checkStringField('custom_addr', false);
            this.checkStringField('custom_port', false);
            this.checkStringField('custom_path', false);
            this.checkStringField('custom_name_eph', false);
            this.checkStringField('custom_name_erp', false);
            this.checkStringField('custom_name_clk', false);            
            
            this.checkStringField('sss_id_list', false);
            this.checkStringField('sss_id_start', false);
            this.checkStringField('sss_id_stop', false);

            this.checkStringField('obs_dir', false, true);
            this.checkCellStringField('obs_name', false);

            this.checkStringField('input_file_ini_path', true);
            
            this.checkStringField('eph_dir', false, true);
            % When the ephemeris file inserted here is not found -> the automatic downloader will dowload the proper file
            this.checkStringField('eph_name', true);           
            this.checkStringField('clk_dir', false, true);
            this.checkStringField('clk_name', true);
            this.checkStringField('crx_dir', false);
            this.checkStringField('dcb_dir', false);
            this.checkStringField('ems_dir', true);

            this.checkStringField('crd_dir', false);
            this.checkStringField('met_dir', false);
            this.checkStringField('ocean_dir', false);

            this.checkStringField('ref_graph_file', true);
            this.checkStringField('geoid_dir', true);
            this.checkStringField('dtm_dir', true);

            this.checkStringField('img_dir', false, true);
            this.checkStringField('img_logo64', false, true);

            this.checkStringField('out_prefix', true);
            this.checkNumericField('run_counter',[0 1e6]);
            
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
                        this.logger.addMessage(sprintf('Checking files for %s receiver %d', obs_type_name, r));
                    end
                    file_name = file_name_all{r};
                    file_count = 0;
                    for f = 1 : numel(file_name)
                        full_path = fnp.checkPath(file_name{f});
                        file_ok = exist(full_path, 'file');
                        file_count = file_count + uint16(logical(file_ok));
                        if go_verbose
                            if file_ok
                                this.logger.addStatusOk(sprintf('%s is present', full_path));
                                status = 1;
                            else
                                this.logger.addError(sprintf('%s does not exist', full_path));
                            end
                        end
                    end
                    if (file_count == 0)
                        status = -1;
                    end
                end                
            end            
        end
        
        function eph_ok = checkNavEphFiles(this, date_start, date_stop)
            % check whether or not all the ephemeris files are available
            eph_ok = true;

            file_name = this.getFullNavEphPath();
                        
            if isempty(file_name)
                eph_ok = false;
            elseif isempty(file_name{1})
                eph_ok = false;
            else
                this.logger.addMarkedMessage('Checking navigational files');
                this.logger.newLine();
                i = 0;
                while (i < numel(file_name) && eph_ok)
                    i = i + 1;
                    eph_ok = exist(file_name{i}, 'file') == 2;
                    if eph_ok
                        this.logger.addStatusOk(sprintf('%s', file_name{i}));
                    else
                        if ~(exist(file_name{i}, 'file') == 7) % if it's not a folder
                            this.logger.addWarning(sprintf('%s does not exist', file_name{i}));
                        else
                            this.logger.addWarning(sprintf('%s it''s a folder, no file name have been declared', file_name{i}));
                        end
                    end
                end
                this.logger.newLine();
            end
        end
        
        function clk_ok = checkNavClkFiles(this)
            % check whether or not all the navigational clock files are available
            
            clk_ok = true;
            file_name = this.getFullNavClkPath();
            
            if isempty(file_name)
                clk_ok = true;
            elseif isempty(file_name{1})
                clk_ok = true;
            else
                this.logger.addMarkedMessage('Checking clock offsets files');
                this.logger.newLine();
                i = 0;
                while (i < numel(file_name) && clk_ok)
                    i = i + 1;
                    clk_ok = exist(file_name{i}, 'file') == 2;
                    if clk_ok
                        this.logger.addStatusOk(sprintf('%s', file_name{i}));
                    else
                        if ~(exist(file_name{i}, 'file') == 7) % if it's not a folder
                            this.logger.addWarning(sprintf('%s does not exist', file_name{i}));
                        else
                            this.logger.addWarning(sprintf('%s it''s a folder, no file name have been declared', file_name{i}));
                        end
                    end
                end
                this.logger.newLine();
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
