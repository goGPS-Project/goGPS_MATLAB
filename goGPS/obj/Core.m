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

%----------------------------------------------------------------------------------------------
%                           goGPS v0.9.1
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
classdef Core < handle
    
    properties (Constant)
    end
    
    properties % Public Access
        cc % Constellation Collector (Enabled Satellite Systems);
    end
       
    methods (Static)
        function obj = Core()
            % Core object creator
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
