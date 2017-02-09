%   CLASS Constellation Collector
% =========================================================================
%
% DESCRIPTION
%   Class to collect and store active Satellite System to be used in the
%   computations
%

%--------------------------------------------------------------------------
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

classdef Constellation_Collector < handle
    properties (Constant, GetAccess = private)
        N_SYS_TOT = 6; % Max number of available satellite systems
        SYS_EXT_NAMES = {'GPS', 'GLONASS', 'Galileo', 'BeiDou', 'QZSS', 'SBAS'}; % full name of the constellation
        SYS_NAMES = {'GPS', 'GLO', 'GAL', 'BDS', 'QZS', 'SBS'}; % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        ID_GPS     = 1 % Id of GPS constellation for goGPS internal use
        ID_GLONASS = 2 % Id of GLONASS constellation for goGPS internal use
        ID_GALILEO = 3 % Id of Galileo constellation for goGPS internal use
        ID_BEIDOU  = 4 % Id of BeiDou constellation for goGPS internal use
        ID_QZSS    = 5 % Id of QZSS constellation for goGPS internal use
        ID_SBAS    = 6 % Id of SBAS constellation for goGPS internal use
        SYS_C      = 'GRECJS'; % Array of constellation char ids GPS = 'G', GLONASS = 'R', ...        
    end
    
    properties (SetAccess = private, GetAccess = public)
        list       % struct of objects (to keep the name of the active constellations)
        sys_name   % rray of active constellations names (as in the list structure)
        num_id     % array of active constellations numeric id
        char_id    % array of active constellations char id
        n_sys      % number of active constellations
        n_sat      % uint8 array of satellite used per active constellation
        n_sat_tot  % uint8 total teoretical number of satellite available for processing
        enabled    % logical array of satellite actually active order: GPS GLO GAL BDS QZS SBS
        
        index      % incremental index of the active satellite system
        prn        % relative id number in the satellite system        
        system     % char id of the constellation per satellite
    end
    
    methods (Access = 'private')
        
        function init(obj, enabled_ss)
            % In it function for the Constellation_Collector Class
            obj.enabled = enabled_ss;
            obj.prn = [];     % relative id number in the satellite system
            obj.system = '';  % char id of the constellation per satellite
            
            obj.n_sat_tot = 0; % counter for number of satellites
            
            if enabled_ss(1) % GPS is active
                obj.list.GPS = GPS_SS(obj.n_sat_tot);
                obj.num_id = [obj.num_id obj.ID_GPS];
                obj.char_id = [obj.char_id obj.list.GPS.char_id];
                obj.system = [obj.system char(ones(1, obj.list.GPS.n_sat) * obj.list.GPS.char_id)];
                obj.prn = [obj.prn; obj.list.GPS.prn];
                obj.n_sat = [obj.n_sat obj.list.GPS.n_sat];
                obj.n_sat_tot = obj.n_sat_tot + obj.list.GPS.n_sat;
            end
            if enabled_ss(2) % GLONASS is active
                obj.list.GLO = GLONASS_SS(obj.n_sat_tot);
                obj.num_id = [obj.num_id obj.ID_GLONASS];
                obj.char_id = [obj.char_id obj.list.GLO.char_id];
                obj.system = [obj.system char(ones(1, obj.list.GLO.n_sat) * obj.list.GLO.char_id)];
                obj.prn = [obj.prn; obj.list.GLO.prn];
                obj.n_sat = [obj.n_sat obj.list.GLO.n_sat];
                obj.n_sat_tot = obj.n_sat_tot + obj.list.GLO.n_sat;
            end
            if enabled_ss(3) % Galileo is active
                obj.list.GAL = Galileo_SS(obj.n_sat_tot);
                obj.num_id = [obj.num_id obj.ID_GALILEO];
                obj.char_id = [obj.char_id obj.list.GAL.char_id];
                obj.system = [obj.system char(ones(1, obj.list.GAL.n_sat) * obj.list.GAL.char_id)];
                obj.prn = [obj.prn; obj.list.GAL.prn];
                obj.n_sat = [obj.n_sat obj.list.GAL.n_sat];
                obj.n_sat_tot = obj.n_sat_tot + obj.list.GAL.n_sat;
            end
            if enabled_ss(4) % BeiDou is active
                obj.list.BDS = BeiDou_SS(obj.n_sat_tot);
                obj.num_id = [obj.num_id obj.ID_BEIDOU];
                obj.char_id = [obj.char_id obj.list.BDS.char_id];
                obj.system = [obj.system char(ones(1, obj.list.BDS.n_sat) * obj.list.BDS.char_id)];
                obj.prn = [obj.prn; obj.list.BDS.prn];
                obj.n_sat = [obj.n_sat obj.list.BDS.n_sat];
                obj.n_sat_tot = obj.n_sat_tot + obj.list.BDS.n_sat;
            end
            if enabled_ss(5) % QZSS is active
                obj.list.QZS = QZSS_SS(obj.n_sat_tot);
                obj.num_id = [obj.num_id obj.ID_QZSS];
                obj.char_id = [obj.char_id obj.list.QZS.char_id];
                obj.system = [obj.system char(ones(1, obj.list.QZS.n_sat) * obj.list.QZS.char_id)];
                obj.prn = [obj.prn; obj.list.QZS.prn];
                obj.n_sat = [obj.n_sat obj.list.QZS.n_sat];
                obj.n_sat_tot = obj.n_sat_tot + obj.list.QZS.n_sat;
            end
            if enabled_ss(6) % SBAS is active (not yet implemented) 
                obj.list.SBAS.n_sat = 0;    % (not yet implemented)
            end        
            
            obj.index = (1 : obj.n_sat_tot)';   % incremental index of the active satellite system
            obj.n_sys = numel(obj.list);
            obj.sys_name = obj.SYS_NAMES(obj.num_id);
        end
        
    end
    
    methods
        function obj = Constellation_Collector(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            % Constructor - parameters: [GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]
            % SYNTAX:
            %   cc = Constellation_Collector('GRECJS'); % use the array of constellation char ids
            %   cc = Constellation_Collector(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
            %   cc = Constellation_Collector([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
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
            %   object handler
            %
            % DESCRIPTION:
            %   Multi-constellation set-up.
            
            % Manually manage overloading
            switch nargin
                case 1
                    if (ischar(GPS_flag))
                        enabled_ss = false(1,obj.N_SYS_TOT);
                        [~, ids] = intersect('GRECJS',GPS_flag);
                        enabled_ss(ids) = true;
                    else
                        enabled_ss = logical(GPS_flag);
                    end
                case 5,  enabled_ss = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, 0]);
                case 6,  enabled_ss = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
                otherwise, error(['Initialization of Constellation_Collector failed: ' 10 '   invalid number of parameters in the constructor call']);
            end
            
            % check the size of the array enabled
            if (numel(enabled_ss) < obj.N_SYS_TOT)
                tmp = false(obj.N_SYS_TOT, 1);
                tmp(1:numel(enabled_ss)) = enabled_ss;
                enabled_ss = tmp;
                clear tmp;
            else
                enabled_ss = enabled_ss(1:obj.N_SYS_TOT);
            end
            
            obj.init(enabled_ss);
        end    
    end

    methods
        function str_cell = toString(obj, str_cell)
            % Display the satellite system in use
            if (nargin == 1)
                str_cell = {};
            end
            [~, ids] = intersect('GRECJS', obj.char_id);
            toString = @(var) regexprep(evalc(['disp(var)']), '''', '');
            str_cell{numel(cell_str) + 1} = ['Constellation in use: ' toString(sort(ids))];
        end
        
        function str_cell = toIniString(obj, str_cell)            
            % Conversion to string of the minimal information needed to reconstruct the obj            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniString('constellations_in_use', obj.char_id, str_cell);
            str_cell = Ini_Manager.toIniString('index', obj.index, str_cell);
            str_cell = Ini_Manager.toIniString('prn', obj.prn, str_cell);
            str_cell = Ini_Manager.toIniString('system', obj.system, str_cell);
        end
        
        
    end
end
