%   CLASS Global_Configuration
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage the parameters of the execution of goGPS
%   This singleton class collects multiple objects containing various
%   parameters, state contains the parameter in use, while the other
%   properties of the class are needed during the execution of goGPS
%
% EXAMPLE
%   settings = Global_Configuration.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Global_Configuration

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
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

classdef Global_Configuration < Settings_Interface

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
        geoid = struct('file', [], 'grid', 0, 'cellsize', 0, 'Xll', 0, 'Yll', 0, 'ncols', 0, 'nrows', 0); % parameters of the reference geoid

        reference = struct('path' , [], 'adj_mat', []);  % reference path for constrained solution, and adjacency matrix

        local_storage = '';
        
        is_advanced = true;
    end

    properties (GetAccess = private, SetAccess = private) % Public Access
        state = Main_Settings();        % Processing settings
    end

    % =========================================================================
    %  INIT
    % =========================================================================
    methods
        % Guard the constructor against external invocation.  We only want
        % to allow a single instance of this class.  See description in
        % Singleton superclass.
        function this = Global_Configuration()
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
    
    methods
        function state = getCurrentSettings(this, ini_settings_file)
            % Get the persistent settings
            %
            % SYNTAX
            %   state = getCurrentSettings(this, ini_settings_file)
            if nargin == 2
                this.state.importIniFile(ini_settings_file);
            end
            % Return the handler to the object containing the current settings
            state = handle(this.state);
        end
        
        function setCurrentSettings(this, cs)
            % Set the persistent settings
            %
            % SYNTAX
            %   this.setCurrentSettings(ini_settings_file)
            this.state = cs;
        end
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods % Public Access
        function import(this, settings)
            % This function try to import settings from another setting object
            if isprop(settings, 'ps')
                this.state.import(settings.state);
            else
                try
                    this.state.import(settings);
                catch ex
                    this.log.addWarning(['Global_Configuration.import failed to import settings (invalid input settings) ', ex.message()]);
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
            str = this.state.toString(str);
        end

        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the obj
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = this.state.export(str_cell);
        end
    end

    % =========================================================================
    %  GOGPS INIT FUNCTIONS
    % =========================================================================
    methods (Access = public)
        function initConfiguration(this)
            % Load all the files necessary to the functioning of a goGPS session
            % SYNTAX:   this.initProcessing()
            
            % Load external resources and update
            fnp = File_Name_Processor();
            out_dir = fnp.checkPath(this.state.getOutDir());
            if ~exist(out_dir, 'dir')
                this.log.newLine();
                this.log.addWarning(sprintf('The out folder does not exists\n Creating %s', out_dir));
                mkdir(out_dir);
            end

            this.log.addMessage(this.log.indent(this.state.cc.toString, 5));

            this.initGeoid();
        end

        function initGeoid(this, geoid)
            % load external geoid
            %
            % SYNTAX
            %   this.initGeoid(); -> load from file
            %   this.initGeoid(geoid); -> import from obj
            if nargin == 2
                this.geoid = geoid;
            else
                try
                    geoid_file = this.state.getGeoidFile();
                    if ~exist(this.state.getGeoidDir, 'file')
                        this.state.geoid_dir = File_Name_Processor.getFullDirPath('../data/reference/geoid', pwd);
                        geoid_file = this.state.getGeoidFile();
                    end
                    if ~strcmp(this.geoid.file, geoid_file)
                        g = load(geoid_file);
                        fn = fieldnames(g);
                        % geoid grid and parameters
                        this.geoid.file = geoid_file;
                        this.geoid.grid = g.(fn{1});
                        this.geoid.cellsize = 360 / size(this.geoid.grid, 2);
                        this.geoid.Xll = -180 + this.geoid.cellsize / 2;
                        this.geoid.Yll = -90 + this.geoid.cellsize / 2;
                        this.geoid.ncols = size(this.geoid.grid, 2);
                        this.geoid.nrows = size(this.geoid.grid, 1);
                        clear g
                    end
                catch
                    this.log.addWarning('Reference geoid not found');
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
        
        function [is_adv] = isAdvanced(this)
            % Get the status of usage (normal/advanced)
            %
            % SYNTAX:
            %   this.isAdvanced()
            is_adv = this.is_advanced;
        end
    end

    % =========================================================================
    %  ADDITIONAL SETTERS
    % =========================================================================
    methods
        function setAdvanced(this, mode)
            % Set the status of usage (normal/advanced)
            %
            % SYNTAX:
            %   this.setAdvanced(<mode>)
            if nargin == 1
                mode = true;
            end
            this.is_advanced = mode;
        end        
    end
    
    % =========================================================================
    %  TEST
    % =========================================================================
    methods (Static, Access = 'public')
        function test()
            % test the class
            s = Global_Configuration.getInstance();
            s.testInterfaceRoutines();
        end
    end
end
