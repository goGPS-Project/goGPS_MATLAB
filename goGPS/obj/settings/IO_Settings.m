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
    end
        
    methods
        function obj = IO_Settings()
            % Creator of IO_settings - verbosity level (true/false) can be set or ini file
        end
    end
        
    methods
        function import(obj, settings)
            % This function import IO (only) settings from another setting object
            if isa(settings, 'Ini_Manager')
                obj.dtm_dir    = settings.getData('dtm_dir');
            else
                obj.dtm_dir    = settings.dtm_dir;
            end
        end
        
        function str = toString(obj, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end            
            str = [str '---- IO SETTINGS ---------------------------------------------------------' 10 10];
            str = [str sprintf(' Directory of DTM folders:                         %s\n\n', obj.dtm_dir)];
        end
        
        function str_cell = export(obj, str_cell)            
            % Conversion to string ini format of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringSection('IO', str_cell);
            str_cell = Ini_Manager.toIniStringComment('Directory of DTM folders', str_cell);
            str_cell = Ini_Manager.toIniString('dtm_dir', obj.dtm_dir, str_cell);
        end
        
    end    
end
