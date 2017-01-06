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
classdef IO_Settings < Settings_Interface
    
    properties (Constant, GetAccess = private)
        DEFAULT_DIR_IN = ['..' filesep 'data' filesep 'data_goGPS'];
        DEFAULT_DIR_OUT = ['..' filesep 'data' filesep 'data_goGPS'];
    end
    
    properties (SetAccess = private, GetAccess = private)        
        logger; % Handler to the logger object

        %------------------------------------------------------------------
        % RECEIVERS
        %------------------------------------------------------------------

        % Observation files of the Receivers
        %receiver_rinex_file = Receiver_File;        
        
        %------------------------------------------------------------------
        % DTM (SET PATH AND LOAD PARAMETER FILES)
        %------------------------------------------------------------------

        % Path to DTM folder containing DTM files
        dtm_dir = ['..' filesep 'data' filesep 'dtm'];
    end
    
    properties (SetAccess = private, GetAccess = public)
        %------------------------------------------------------------------
        % DTM (SET PATH AND LOAD PARAMETER FILES)
        %------------------------------------------------------------------
        
        % Parameters common to all DTM tiles
        dtm_tile_header = struct('nrows', 0, 'ncols', 0, 'cellsize', 0, 'nodata', 0);
        
        % Parameters used to georeference every DTM tile
        dtm_tile_georef = zeros(1,1,4);
    end
        
    methods
        function obj = IO_Settings()
            % Creator of IO_settings - verbosity level (true/false) can be set
            obj.logger = Logger.getInstance();
            obj.init_dtm();
        end
    end
    
    methods 
        
        function init_dtm(obj, dtm_dir)
            % Try to load default DTM values
            if nargin == 1
                dtm_dir = obj.dtm_dir;
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
    
    methods
        function copyFrom(obj, io_settings)
            % This function import IO (only) settings from another setting object
            obj.dtm_dir    = io_settings.dtm_dir;
            obj.init_dtm();
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- IO SETTINGS ---------------------------------------------------------' 10 10];
            str = [str sprintf(' directory of DTM:                        %s\n\n', obj.dtm_dir)];
        end
        
        function str_cell = toIniString(obj, str_cell)            
            % Conversion to string of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringSection('IO', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', obj.dtm_dir, str_cell);
        end
        
    end    
end
