%   CLASS Meteo_Network
% =========================================================================
%
% DESCRIPTION
%   Class to store and manage a network of meteo station
%
% EXAMPLE
%   settings = Meteo_Network();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Meteo_Data
%
% REFERENCE

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0RC1
%
%--------------------------------------------------------------------------
%  Copyright (C) 2021 (GReD srl) Andrea Gatti, Giulio Tagliaferro, Eugenio Realini
%  Written by:        Andrea Gatti, Giulio Tagliaferro ...
%  Contributors:      Andrea Gatti, Giulio Tagliaferro ...
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
%----------------------------------------------------------------------------------------------
classdef Meteo_Network < handle
    properties
        mds % list of meteo data
        log % handfle of logger
    end
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static)
        % Concrete implementation.  See Singleton superclass.
        function this = Meteo_Network()
            this.log = Core.getLogger();
        end
    end
    
    %% METHODS UI
    % ==================================================================================================================================================
    methods (Static, Access = public)    
        function ok_go = openGUI()
            ok_go = gui_goGPS;
        end
    end
    
    methods (Access = public)
        function initSession(this, date_start, date_stop)
            % Load all the metereological data
            %
            % SYNTAX
            % this.initSession(data_start, data_stop)
 
            this.mds = [];
            % load all meteo file present in current settings
            state = Core.getCurrentSettings();
            fnames = state.getMetFileName(date_start, date_stop);
            this.initFromFileList(fnames)
        end
        
        function initFromFileList(this, file_name_list)
            % Load metereological data directly from file list
            %
            % INPUT
            %   file_name_list: two level cells, first level station second level days
            %
            % SYNTAX
            %   this.initFromFileList(file_name_list)
            
            n_stations = numel(file_name_list);
            for  i = 1 : n_stations
                n_sss = numel(file_name_list{i});
                for  s = 1 : n_sss
                    if ~isempty(File_Name_Processor.getFileName(file_name_list{i}{s}))
                        if ~exist(file_name_list{i}{s}, 'file')
                            this.log.addWarning(sprintf('Skipping %s - file not found', file_name_list{i}{s}));
                        else
                            md = Meteo_Data(file_name_list{i}{s});
                            if md.isValid()
                                % md.setMaxBound(0);
                                if length(this.mds) < i
                                    this.mds = [this.mds; md];
                                else
                                    this.mds(i).inject(md);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function md = getVMS(this, name, xyz, time)
            % Get Virtual Meteo Station
            %
            % INPUT
            %   name    name of the new Meteo_Data virtual station
            %   xyz     coordinates of the new meteo station
            %   time    time of interpolation
            %
            %   
            % OUTPUT
            %   md      virtual Meteo_Data station generated at xyz coordinates
            %   
            % SYNTAX
            %   md = this.getVMS(name, xyz, time)
            %
            % EXAMPLE
            %   [x, y, z, amsl] = station(1).getLocation();
            %   md1 = this.getVMS('test', [x y z], station(1).getObsTime)
            md = Meteo_Data.getVMS( name, xyz, time, this.mds);
        end
    end
    
end
