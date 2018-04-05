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
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.6.0 alpha 2 - nightly
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2018 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:     Giulio Tagliaferro, ...
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
    end
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods (Static, Access = private)
        % Concrete implementation.  See Singleton superclass.
        function this = Meteo_Network()
            
        end
    end
    
    %% METHODS UI
    % ==================================================================================================================================================
    methods (Static, Access = public)
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_meteo_network__

            if isempty(unique_instance_meteo_network__)
                this = Meteo_Network();
                unique_instance_meteo_network__ = this;
            else
                this = unique_instance_meteo_network__;
            end
            
        end

        function ok_go = openGUI()
            ok_go = gui_goGPS;
        end
    end
    
    methods (Access = public)
        function initSession(this, date_start, date_stop)
            this.mds = [];
            % load all meteo file present in current settings
            state = Global_Configuration.getCurrentSettings();
            fnames = state.getMetFileName(date_start, date_stop);
            n_met_data = numel(fnames);
            for  i=1:n_met_data
                
                md = Meteo_Data(fnames{i});
                if md.isValid()
                    %md.setMaxBound(0);
                    this.mds = [this.mds; md];
                
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