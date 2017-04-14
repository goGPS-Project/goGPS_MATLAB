%   CLASS Core
% =========================================================================
%
% DESCRIPTION
%   Collector of settings to manage a useful parameters of goGPS
%   This singleton class collects multiple objects containing various
%   parameters
%
% EXAMPLE
%   settings = goGNSS();
%
% FOR A LIST OF CONSTANTs and METHODS use doc goGNSS


%--------------------------------------------------------------------------
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

classdef Core < handle
    
    properties (Constant)
        GO_GPS_VERSION = '0.5.1 beta';
    end
    
    properties % Public Access
        logger;        
    end
    
    methods (Static)
        function this = Core()
            % Core object creator
            this.logger = Logger.getInstance();
        end
        
        % Concrete implementation.  See Singleton superclass.
        function this = getInstance()
            % Get the persistent instance of the class
            persistent unique_instance_core__
            unique_instance_core__ = [];
            
            if isempty(unique_instance_core__)
                this = Core();
                unique_instance_core__ = this;
            else
                this = unique_instance_core__;
            end
        end
        
        
        function showTextHeader()
            this.logger = Logger.getInstance();
            if this.logger.getColorMode()
                cprintf([241 160 38]/255,'               ___ ___ ___\n     __ _ ___ / __| _ | __|\n    / _` / _ \\ (_ |  _|__ \\\n    \\__, \\___/\\___|_| |___/\n    |___/                    '); cprintf('text','v '); cprintf('text', Core.GO_GPS_VERSION); fprintf('\n');
                fprintf('\n--------------------------------------------------------------------------\n\n');
            else                
                fprintf('--------------------------------------------------------------------------\n');
                fprintf('               ___ ___ ___\n');
                fprintf('     __ _ ___ / __| _ | __|\n');
                fprintf('    / _` / _ \ (_ |  _|__ \\n');
                fprintf('    \__, \___/\___|_| |___/\n');
                fprintf('    |___/                    v %s\n', Core.GO_GPS_VERSION);
                fprintf('\n');
                fprintf('--------------------------------------------------------------------------\n\n');
            end
        end
    end
    % =========================================================================
    %  CONSTELLATION MANAGEMENT
    % =========================================================================
    
    methods % Public Access
        function [ cc ] = initConstellation(obj, GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            % Multi-constellation set-up.
            %
            % SYNTAX:
            %   cc = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
            %   cc = initConstellation([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
            %
            % INPUT:
            %   single logical array whose elements are:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   the results is stored within the object referenced by "cc"
            switch nargin
                case 2,  enabled_ss = GPS_flag;
                case 6,  enabled_ss = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, 0]);
                case 7,  enabled_ss = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
                otherwise, error(['Initialization of Constellation_Collector failed: ' 10 '   invalid number of parameters in the constructor call']);
            end
            
            obj.cc = Costellation_Collector(enabled_ss);
            
            cc = obj.cc;
        end
        
        function [cc] = getConstellations(obj)
            % Get the object containing the actual status of the enabled constallation
            cc = obj.cc;
        end
    end
    
    methods % Public Access (Legacy support)
    end
    
end
