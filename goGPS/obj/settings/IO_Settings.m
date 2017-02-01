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
%    |___/                    v 0.9.1
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
    
    properties (Constant, GetAccess = private)
        DEFAULT_DIR_IN = ['..' filesep 'data' filesep];
        DEFAULT_DIR_OUT = ['..' filesep 'data' filesep];  
    end
    
    properties (SetAccess = private, GetAccess = private)
        %------------------------------------------------------------------
        % RECEIVERS
        %------------------------------------------------------------------

        % Observation files of the Receivers
        %receiver_rinex_file = Receiver_File;        

        %------------------------------------------------------------------
        % SATELLITES
        %------------------------------------------------------------------
        
        % Path to Navigational files folder
        nav_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'NAV' filesep];
        
        % Path to Clock Offset files folder
        clk_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'CLK' filesep];
        
        % Path to CRX folder containing files of Satellites problems
        crx_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'CRX' filesep];

        % Path to DCB folder containing files of Differential Code Biases
        dcb_dir = [IO_Settings.DEFAULT_DIR_IN 'SAT' filesep 'DCB' filesep];

        %------------------------------------------------------------------
        % REFERENCE
        %------------------------------------------------------------------
        
        % Path to Geoid folder containing the geoid to be used for the computation of hortometric heighs
        geoid_dir = [IO_Settings.DEFAULT_DIR_IN 'geoid' filesep];
        
        %------------------------------------------------------------------
        % DTM (SET PATH AND LOAD PARAMETER FILES)
        %------------------------------------------------------------------

        % Path to DTM folder containing DTM files
        dtm_dir = [IO_Settings.DEFAULT_DIR_IN 'dtm' filesep];
        
        %------------------------------------------------------------------
        % UI IMAGES
        %------------------------------------------------------------------

        % Path to images used by the interface
        img_dir = [IO_Settings.DEFAULT_DIR_IN 'img' filesep];
        
    end
    
    properties (SetAccess = private, GetAccess = public)
    end
        
    methods
        function this = IO_Settings()
            % Creator of IO_settings - verbosity level (true/false) can be set or ini file
        end
    end
        
    methods
        function import(this, settings)
            % This function import IO (only) settings from another setting object
            if isa(settings, 'Ini_Manager')
                this.nav_dir    = settings.getData('nav_dir');
                this.clk_dir    = settings.getData('clk_dir');
                this.crx_dir    = settings.getData('crx_dir');
                this.dcb_dir    = settings.getData('dcb_dir');
                this.geoid_dir  = settings.getData('geoid_dir');
                this.dtm_dir    = settings.getData('dtm_dir');
                this.img_dir    = settings.getData('img_dir');
            else
                this.nav_dir    = settings.nav_dir;
                this.clk_dir    = settings.clk_dir;
                this.crx_dir    = settings.crx_dir;
                this.dcb_dir    = settings.dcb_dir;
                this.geoid_dir  = settings.geoid_dir;
                this.dtm_dir    = settings.dtm_dir;
                this.img_dir    = settings.img_dir;
            end
        end
        
        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- IO SETTINGS ---------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of Navigational Files:                  %s\n', this.nav_dir)];
            str = [str sprintf(' Directory of Satellite clock offsets:             %s\n', this.clk_dir)];
            str = [str sprintf(' Directory of CRX (satellite problems):            %s\n', this.crx_dir)];
            str = [str sprintf(' Directory of DCB (Differential Code Biases):      %s\n\n', this.dcb_dir)];
            str = [str sprintf(' Directory of Geoid models:                        %s\n\n', this.geoid_dir)];
            str = [str sprintf(' Directory of DTM data:                            %s\n\n', this.dtm_dir)];
            str = [str sprintf(' Directory of images for UI:                       %s\n\n', this.img_dir)];
        end
        
        function str_cell = export(this, str_cell)            
            % Conversion to string ini format of the minimal information needed to reconstruct the this            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringSection('IO', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Navigational files', str_cell);
            str_cell = Ini_Manager.toIniString('nav_dir', this.nav_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of clock offset files', str_cell);
            str_cell = Ini_Manager.toIniString('crx_dir', this.clk_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of CRX files (containing satellite problems)', str_cell);
            str_cell = Ini_Manager.toIniString('crx_dir', this.crx_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DCB files (Differential Code Biases)', str_cell);
            str_cell = Ini_Manager.toIniString('dcb_dir', this.dcb_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of Geoid files', str_cell);
            str_cell = Ini_Manager.toIniString('geoid_dir', this.geoid_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DTM data', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', this.dtm_dir, str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of images for UI', str_cell);
            str_cell = Ini_Manager.toIniString('img_dir', this.img_dir, str_cell);
        end
        
    end    
end
