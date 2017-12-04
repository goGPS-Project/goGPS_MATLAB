%   CLASS Go_State
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage the parameters of the execution of goGPS
%   This singleton class collects multiple objects containing various
%   parameters, cur_settings contains the parameter in use, while the other
%   properties of the class are needed during the execution of goGPS
%
% EXAMPLE
%   settings = Go_State.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Go_State

%--------------------------------------------------------------------------
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

classdef Go_State < Settings_Interface

    properties (Constant)
        V_LIGHT = 299792458;                % Velocity of light in the void [m/s]

        % Values as defined by standards

        PI_ORBIT = 3.1415926535898;         % pi as from standards
        CIRCLE_RAD = 6.2831853071796;       % Circle as from standards

        % Standard atmosphere parameters (struct PRES, STD_TEMP. STD_HUMI) - Berg, 1948
        ATM = struct('PRES', 1013.25, ...      % pressure [mbar]
                     'STD_TEMP', 291.15, ...   % temperature [K] (18? C)
                     'STD_HUMI', 50.0);        % humidity [%]
    end

    properties (GetAccess = private, SetAccess = private) % Public Access
        geoid;                                           % parameters of the reference geoid

        reference = struct('path' , [], 'adj_mat', []);  % reference path for constrained solution, and adjacency matrix

        local_storage = '';
    end

    properties % Public Access
        cur_settings = Main_Settings();        % Processing settings
    end

    % =========================================================================
    %  INIT
    % =========================================================================
    methods (Access = private)
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function this = Go_State()
            this.initLogger();

            if ispc()
                home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
                this.local_storage = [home '\AppData\Local\goGPS'];
            else
                home = getenv('HOME');
                if ismac()
                    this.local_storage = [home '/Library/Application Support/goGPS'];
                else
                    this.local_storage = [home '/.goGPS'];
                end
            end
            if ~(exist(this.local_storage, 'dir'))
                mkdir(this.local_storage)
            end
        end
    end

    % =========================================================================
    %  SINGLETON GETTERS
    % =========================================================================
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function this = getInstance(ini_settings_file, force_clean)
            % Get the persistent instance of the class
            persistent unique_instance_settings__
            ini_is_present = false;
            switch nargin
                case 0
                    force_clean = false;
                case 1
                    if ischar(ini_settings_file)
                        ini_is_present = true;
                        force_clean = false;
                    else
                        force_clean = ini_settings_file;
                    end
                case 2
                    ini_is_present = true;
            end

            if force_clean
                log = Logger.getInstance();
                log.addWarning('Cleaning settings of the session');
                clear unique_instance_settings__;
                unique_instance_settings__ = [];
            end

            if isempty(unique_instance_settings__)
                if ini_is_present
                    this = Go_State(ini_settings_file);
                else
                    this = Go_State();
                end
                unique_instance_settings__ = this;
            else
                if ini_is_present
                    this = unique_instance_settings__;
                    this.cur_settings.importIniFile(ini_settings_file);
                else
                    this = unique_instance_settings__;
                    if isempty(this.cur_settings.ext_ini)
                        this.cur_settings.updateExternals();
                    end
                end
            end
        end

        function cur_settings = getCurrentSettings(ini_settings_file)
            % Get the persistent sittings
            if nargin == 0
                this = Go_State.getInstance();
            else
                this = Go_State.getInstance(ini_settings_file);
            end
            % Return the handler to the object containing the current settings
            cur_settings = handle(this.cur_settings);
        end
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods % Public Access
        function import(this, settings)
            % This function try to import settings from another setting object
            if isprop(settings, 'ps')
                this.cur_settings.import(settings.cur_settings);
            else
                try
                    this.cur_settings.import(settings);
                catch ex
                    this.log.addWarning(['Go_State.import failed to import settings (invalid input settings) ', ex.message()]);
                end
            end
        end

        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            str = [str '---- CONSTANTS -----------------------------------------------------------' 10 10];
            str = [str sprintf(' VLIGHT:                                           %g\n', this.V_LIGHT)];
            str = [str sprintf(' PI_ORBIT:                                         %g\n', this.PI_ORBIT)];
            str = [str sprintf(' CIRCLE_RAD:                                       %g\n', this.CIRCLE_RAD)];
            str = [str sprintf(' STANDARD ATMOSPHERE (Berg, 1948):\n  - PRESSURE [mBar]                                %g\n  - TEMPERATURE [K]                                %g\n  - HUMIDITY [%%]                                   %g\n\n', this.ATM.PRES, this.ATM.STD_TEMP, this.ATM.STD_HUMI)];
            str = this.cur_settings.toString(str);
        end

        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = this.cur_settings.export(str_cell);
        end
    end

    % =========================================================================
    %  GOGPS INIT FUNCTIONS
    % =========================================================================
    methods (Access = public)
        function initProcessing(this)
            % Load all the files necessary to the functioning of a goGPS session
            % SYNTAX:   this.initProcessing()
            
            % Load external resources and update
            fnp = File_Name_Processor();
            out_dir = fnp.checkPath(this.cur_settings.getOutDir());
            if ~exist(out_dir, 'dir')
                this.log.newLine();
                this.log.addWarning(sprintf('The out folder does not exists\n Creating %s', out_dir));
                mkdir(out_dir);
            end

            this.log.addMessage(this.log.indent(this.cur_settings.cc.toString, 5));

            this.initRef();
            this.initGeoid();
            
            this.log.addMarkedMessage('Conjuring files!');
            fw = File_Wizard;
            c_mode = this.log.getColorMode();
            this.log.setColorMode(0);
            fw.conjureFiles();
            this.log.setColorMode(c_mode);
        end

        function initGeoid(this)
            % load external geoid
            % SYNTAX: this.initGeoid();
            try
                geoid_file = this.cur_settings.getGeoidFile();
                g = load(geoid_file);
                fn = fieldnames(g);
                % geoid grid and parameters
                this.geoid.grid = g.(fn{1});
                this.geoid.cellsize = 360 / size(this.geoid.grid, 2);
                this.geoid.Xll = -180 + this.geoid.cellsize / 2;
                this.geoid.Yll = -90 + this.geoid.cellsize / 2;
                this.geoid.ncols = size(this.geoid.grid, 2);
                this.geoid.nrows = size(this.geoid.grid, 1);
                clear g
            catch
                this.log.addWarning('Reference geoid not found', 50);
                % geoid unavailable
                this.geoid.grid = 0;
                this.geoid.cellsize = 0;
                this.geoid.Xll = 0;
                this.geoid.Yll = 0;
                this.geoid.ncols = 0;
                this.geoid.nrows = 0;
            end
        end

    end

    methods (Access = private)
        function initRef(this)
            % load external ref_path

            %-------------------------------------------------------------------------------------------
            % REFERENCE PATH LOAD
            %-------------------------------------------------------------------------------------------
            if this.cur_settings.plot_ref_path
                filename_ref = this.cur_settings.getRefFile();
                d = dir(filename_ref);

                if ~isempty(d)
                    load(filename_ref, 'ref_path', 'mat_path');

                    % adjust the reference path according to antenna height
                    [ref_phi, ref_lam, ref_h] = cart2geod(ref_path(:,1),ref_path(:,2),ref_path(:,3)); %#ok<NODEF>
                    ref_h = ref_h + this.cur_settings.antenna_h;
                    orbital_p = this.cur_settings.cc.gps.ORBITAL_P; % Using GPS orbital parameters
                    [ref_X, ref_Y, ref_Z] = geod2cart(ref_phi, ref_lam, ref_h, orbital_p.ELL.A, orbital_p.ELL.F);
                    this.reference.path = [ref_X , ref_Y , ref_Z];
                    this.reference.adj_mat = mat_path;
                else
                    this.reference.path = [];
                    this.reference.adj_mat = [];
                end

            else
                this.reference.path = [];
                this.reference.adj_mat = [];
            end
        end

    end

    % =========================================================================
    %  ADDITIONAL GETTERS
    % =========================================================================
    methods
        function [go_dir] = getLocalStorageDir(this)
            % Get local storage
            go_dir = this.local_storage;
        end

        function [ref_path, mat_path] = getReferencePath(this)
            % Get reference path
            if (nargout == 2)
                ref_path = this.reference.path;
                mat_path = this.reference.adj_mat;
            elseif (nargout == 1)
                ref_path = this.reference;
            end

        end

        function [geoid] = getRefGeoid(this)
            % Get reference path
            if isempty(this.geoid)
                this.initGeoid();
            end
            geoid = this.geoid;
        end
    end

    % =========================================================================
    %  GOGPS EXPORT
    % =========================================================================
    methods
        function varargout = settingsToGo(this, state)
            % function that simulate the interface returns value
            % legacy approach, here the object creates various global variables needed for the processing
            % SYNTAX: varargout = this.settingsToGo(state)
            
            %initGeoid();

            % export settings as they are exported by the GUI
            if nargin == 1
                state = this.cur_settings;
            end

            varargout = cell(42,1);
            varargout{1}  = state.getMode();         % mode
            varargout{2}  = state.constrain;       % constrain
            varargout{3}  = 0;                       % it was mode_data, now goGPS bin files are unsupported, dropping support
            varargout{4}  = state.plot_ref_path;   % plot ref_path
            varargout{5}  = state.flag_rinex_mpos; % get Master position from RINEX
            varargout{6}  = state.plot_master;     % plot master
            varargout{7}  = state.plot_google_earth; % plot Google Earth
            varargout{8}  = state.plot_err_ellipse;  % plot error ellipse
            varargout{9}  = state.flag_ntrip;      % flag NTRIP
            varargout{10} = state.plot_ambiguities;  % plot ambiguities
            varargout{11} = state.plot_skyplot_snr; % plot skyplot
            varargout{12} = state.plot_proc;   % plot while processing
            varargout{13} = state.isVariableKF();   % flag kalman filter mode variable
            varargout{14} = state.stop_go_stop; % stop go stop mode
            varargout{15} = state.cc.isSbsActive(); % use sbas
            varargout{16} = state.flag_iar; % use iar

            file_root_out = checkPath([state.out_dir filesep state.out_prefix '_' num2str(state.run_counter,'%03d') ]);

            if state.isModePP()
                file_name_nav = '';
                file_name_ref = state.getReferencePath();
                file_name_pco = state.getAtxFile();
                file_name_blq = state.getOceanFile();

                if(state.isModeMultiReceiver())
                    [multi_antenna_rf, ~] = state.getGeometry();
                else
                    multi_antenna_rf = [];
                end

            else
                file_name_nav = '';
                file_name_ref = '';
                file_name_pco = '';
                file_name_blq = '';
                multi_antenna_rf = [];
            end

            protocol_idx = state.c_prtc;


            varargout{17} = '';      % prefix of the in binary file no more supported
            varargout{18} = file_root_out;     % prefix of the out prefix
            if state.isModeSEID()
                varargout{19} = state.getReferencePath();
                varargout{20} = state.getTargetPath();
            else
                varargout{19} = state.getTargetPath();
                varargout{20} = state.getMasterPath();
            end
            varargout{21} = file_name_nav;
            varargout{22} = file_name_ref;
            varargout{23} = file_name_pco;
            varargout{24} = file_name_blq;
            varargout{25} = [ state.mpos.X; state.mpos.Y; state.mpos.Z]; % position master station
            varargout{26} = protocol_idx;
            varargout{27} = multi_antenna_rf;

            varargout{28} = state.iono_model;      % iono model
            varargout{29} = state.tropo_model;      % tropo

            % mixed
            varargout{30} = 'default';

            varargout{31} = state.flag_ocean;
            varargout{32} = state.flag_outlier;
            varargout{33} = state.flag_outlier_OLOO;
            varargout{34} = state.flag_tropo;
            varargout{35} = state.flag_tropo_gradient;
            varargout{36} = find(state.cc.getGPS().flag_f);
            varargout{37} = state.isModeSEID();
            varargout{38} = state.p_rate;
            varargout{39} = iif(state.flag_ionofree,'IONO_FREE', 'NONE');
            varargout{40} = state.flag_pre_pro;       % flag pre-processing
            varargout{41} = state.getCrdFile();
            varargout{42} = state.getMetFile();

            global sigmaq0 sigmaq_vE sigmaq_vN sigmaq_vU sigmaq_vel
            global sigmaq_cod1 sigmaq_cod2 sigmaq_codIF sigmaq_ph sigmaq_phIF sigmaq0_N sigmaq_dtm sigmaq0_tropo sigmaq0_tropo_gradient sigmaq_tropo sigmaq_tropo_gradient sigmaq0_rclock sigmaq_rclock
            global min_nsat min_arc cutoff snr_threshold cs_threshold_preprocessing cs_threshold weights snr_a snr_0 snr_1 snr_A order o1 o2 o3
            global h_antenna
            global tile_header tile_georef dtm_dir
            global master_ip master_port ntrip_user ntrip_pw ntrip_mountpoint
            global nmea_init
            global flag_doppler_cs
            global COMportR
            global IAR_method P0 mu flag_auto_mu flag_default_P0
            global SPP_threshold max_code_residual max_phase_residual

            IAR_method = state.iar_mode;
            P0 = state.iar_p0;
            mu = state.iar_mu;
            flag_auto_mu = state.flag_iar_auto_mu;
            flag_default_P0 = state.flag_iar_default_p0;

            SPP_threshold = state.pp_spp_thr;
            max_code_residual = state.pp_max_code_err_thr;
            max_phase_residual = state.pp_max_phase_err_thr;


            if state.c_n_receivers >= 1
                COMportR{1,1} = state.c_com_addr{1};
                if state.c_n_receivers >= 2
                    COMportR{2,1} = state.c_com_addr{1};
                    if state.c_n_receivers >= 3
                        COMportR{3,1} = state.c_com_addr{1};
                        if state.c_n_receivers >= 4
                            COMportR{4,1} = state.c_com_addr{1};
                        end
                    end
                end
            end

            cutoff = state.cut_off;
            snr_threshold = state.snr_thr;
            flag_doppler_cs = state.flag_doppler;
            sigmaq0 = state.sigma0_k_pos ^ 2;
            sigmaq_vE = state.std_k_ENU.E ^ 2;
            sigmaq_vN = state.std_k_ENU.N ^ 2;
            sigmaq_vU = state.std_k_ENU.U ^ 2;
            sigmaq_vel = state.std_k_vel_mod ^ 2;
            sigmaq_cod1 = state.std_code ^  2;
            sigmaq_cod2 = 0.16;                 % <-- to be changed in the future
            sigmaq_codIF = 1.2 ^ 2;             % <-- to be changed in the future

            sigmaq_ph = state.std_phase ^ 2;
            sigmaq_phIF = iif(state.std_phase == 1e30, 1e30, 0.009^2); % <-- to be changed in the future
            sigmaq0_N = 1000;

            sigmaq_dtm = state.std_dtm ^ 2;
            sigmaq0_tropo = state.sigma0_tropo ^ 2;
            sigmaq0_tropo_gradient = state.sigma0_tropo_gradient ^ 2;
            sigmaq_tropo = state.std_tropo ^ 2;               %(0.005/sqrt(120))^2
            sigmaq_tropo_gradient = state.std_tropo_gradient ^ 2;    %(0.001/sqrt(120))^2
            sigmaq0_rclock = 2e-17;             % <-- to be changed in the future
            sigmaq_rclock = 1e3;                % <-- to be changed in the future

            min_nsat = state.min_n_sat;
            min_arc = state.min_arc;

            cs_threshold = state.cs_thr;
            cs_threshold_preprocessing = state.cs_thr_pre_pro;

            weights = state.w_mode;
            snr_a = state.w_snr.a;
            snr_0 = state.w_snr.zero;
            snr_1 = state.w_snr.one;
            snr_A = state.w_snr.A;

            global amb_restart_method
            amb_restart_method = state.iar_restart_mode;
            order = state.kf_mode + 1;

            o1 = order;
            o2 = order*2;
            o3 = order*3;

            h_antenna = state.antenna_h;

            dtm_dir = state.dtm_dir;
            try
                load([dtm_dir '/tiles/tile_header'], 'tile_header');
                load([dtm_dir '/tiles/tile_georef'], 'tile_georef');
            catch e
                tile_header.nrows = 0;
                tile_header.ncols = 0;
                tile_header.cellsize = 0;
                tile_header.nodata = 0;
                tile_georef = zeros(1,1,4);
            end
            master_ip = state.ntrip.ip_addr;
            master_port = iif(ischar(state.ntrip.port), str2num(state.ntrip.port), state.ntrip.port);
            ntrip_user = state.ntrip.username;
            ntrip_pw = state.ntrip.password;
            ntrip_mountpoint = state.ntrip.mountpoint;
            phiApp = state.ntrip.approx_position.lat;
            lamApp = state.ntrip.approx_position.lon;
            hApp = state.ntrip.approx_position.h;
            [XApp,YApp,ZApp] = geod2cart (phiApp*pi/180, lamApp*pi/180, hApp, 6378137, 1/298.257222101);
            if ~isnan(XApp) && ~isnan(YApp) && ~isnan(ZApp)
                nmea_init = NMEA_GGA_gen([XApp YApp ZApp],10);
            else
                nmea_init = '';
            end
        end
    end

    % =========================================================================
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()
            % test the class
            s = Go_State.getInstance();
            s.testInterfaceRoutines();
        end
    end
end
