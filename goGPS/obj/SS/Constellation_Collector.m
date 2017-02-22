%   CLASS Constellation Collector
% =========================================================================
%
% DESCRIPTION
%   Class to collect and store active Satellite System to be used in the
%   computations
%
% EXAMPLE
%   cc = Constellation_Collector('GRECJS');
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

classdef Constellation_Collector < Settings_Interface
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
        sys_name   % array of active constellations names (as in the list structure)
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
        function init(this, enabled_ss)
            % In it function for the Constellation_Collector Class
            if (sum(enabled_ss) > 0)
                this.enabled = enabled_ss;
                this.prn = [];     % relative id number in the satellite system
                this.system = '';  % char id of the constellation per satellite
                
                this.n_sat_tot = 0; % counter for number of satellites
                
                if enabled_ss(1) % GPS is active
                    this.list.GPS = GPS_SS(this.n_sat_tot);
                    this.num_id = [this.num_id this.ID_GPS];
                    this.char_id = [this.char_id this.list.GPS.char_id];
                    this.system = [this.system char(ones(1, this.list.GPS.n_sat) * this.list.GPS.char_id)];
                    this.prn = [this.prn; this.list.GPS.prn];
                    this.n_sat = [this.n_sat this.list.GPS.n_sat];
                    this.n_sat_tot = this.n_sat_tot + this.list.GPS.n_sat;
                end
                if enabled_ss(2) % GLONASS is active
                    this.list.GLO = GLONASS_SS(this.n_sat_tot);
                    this.num_id = [this.num_id this.ID_GLONASS];
                    this.char_id = [this.char_id this.list.GLO.char_id];
                    this.system = [this.system char(ones(1, this.list.GLO.n_sat) * this.list.GLO.char_id)];
                    this.prn = [this.prn; this.list.GLO.prn];
                    this.n_sat = [this.n_sat this.list.GLO.n_sat];
                    this.n_sat_tot = this.n_sat_tot + this.list.GLO.n_sat;
                end
                if enabled_ss(3) % Galileo is active
                    this.list.GAL = Galileo_SS(this.n_sat_tot);
                    this.num_id = [this.num_id this.ID_GALILEO];
                    this.char_id = [this.char_id this.list.GAL.char_id];
                    this.system = [this.system char(ones(1, this.list.GAL.n_sat) * this.list.GAL.char_id)];
                    this.prn = [this.prn; this.list.GAL.prn];
                    this.n_sat = [this.n_sat this.list.GAL.n_sat];
                    this.n_sat_tot = this.n_sat_tot + this.list.GAL.n_sat;
                end
                if enabled_ss(4) % BeiDou is active
                    this.list.BDS = BeiDou_SS(this.n_sat_tot);
                    this.num_id = [this.num_id this.ID_BEIDOU];
                    this.char_id = [this.char_id this.list.BDS.char_id];
                    this.system = [this.system char(ones(1, this.list.BDS.n_sat) * this.list.BDS.char_id)];
                    this.prn = [this.prn; this.list.BDS.prn];
                    this.n_sat = [this.n_sat this.list.BDS.n_sat];
                    this.n_sat_tot = this.n_sat_tot + this.list.BDS.n_sat;
                end
                if enabled_ss(5) % QZSS is active
                    this.list.QZS = QZSS_SS(this.n_sat_tot);
                    this.num_id = [this.num_id this.ID_QZSS];
                    this.char_id = [this.char_id this.list.QZS.char_id];
                    this.system = [this.system char(ones(1, this.list.QZS.n_sat) * this.list.QZS.char_id)];
                    this.prn = [this.prn; this.list.QZS.prn];
                    this.n_sat = [this.n_sat this.list.QZS.n_sat];
                    this.n_sat_tot = this.n_sat_tot + this.list.QZS.n_sat;
                end
                if enabled_ss(6) % SBAS is active (not yet implemented)
                    this.list.SBAS.n_sat = 0;   % (not yet implemented)
                    this.num_id = [this.num_id this.ID_SBAS];
                end
                
                this.index = (1 : this.n_sat_tot)';   % incremental index of the active satellite system
                this.n_sys = numel(this.num_id);
                this.sys_name = this.SYS_NAMES(this.num_id);
            else
                this.logger.addError('No satellite system selected -> Enabling GPS');
                ss = false(this.N_SYS_TOT,1); ss(1) = true;
                this.init(ss);
            end
        end
    end
    
    methods
        function this = Constellation_Collector(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
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
                        enabled_ss = false(1,this.N_SYS_TOT);
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
            if (numel(enabled_ss) < this.N_SYS_TOT)
                tmp = false(this.N_SYS_TOT, 1);
                tmp(1:numel(enabled_ss)) = enabled_ss;
                enabled_ss = tmp;
                clear tmp;
            else
                enabled_ss = enabled_ss(1:this.N_SYS_TOT);
            end
            
            this.init(enabled_ss);
        end    
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================    
    methods (Access = 'public')
        function import(this, settings)
            % This function import processing settings from another setting object or ini file
            enabled_ss = false(1,this.N_SYS_TOT);
            if isa(settings, 'Ini_Manager')
                [~, ids] = intersect(this.SYS_C,settings.getData('active_constellation_ch'));
            else
                [~, ids] = intersect(this.SYS_C,settings.char_id);
            end
            enabled_ss(ids) = true;
            this.init(enabled_ss);
        end
        
        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            [~, ids] = intersect(this.SYS_C, this.char_id);
            toString = @(var) regexprep(regexprep(evalc(['disp(var)']), '''   ', ','), '''', '');
            str = [str ' Constellation in use: ' toString(this.SYS_EXT_NAMES(sort(ids)))];
        end
        
        function str_cell = export(this, str_cell)
            % Conversion to string of the minimal information needed to reconstruct the this            
            if (nargin == 1)
                str_cell = {};
            end
            str_cell = Ini_Manager.toIniStringComment('Active constallations for the processing', str_cell);
            str_cell = Ini_Manager.toIniStringComment('insert the constellations as array of char:', str_cell);
            str_cell = Ini_Manager.toIniStringComment(' - "G" GPS', str_cell);
            str_cell = Ini_Manager.toIniStringComment(' - "R" GLONASS', str_cell);
            str_cell = Ini_Manager.toIniStringComment(' - "E" Galileo', str_cell);
            str_cell = Ini_Manager.toIniStringComment(' - "C" BeiDou', str_cell);
            str_cell = Ini_Manager.toIniStringComment(' - "J" QZSS', str_cell);
            str_cell = Ini_Manager.toIniStringComment(' - "S" SBAS (not yet available)', str_cell);
            str_cell = Ini_Manager.toIniString('active_constellation_ch', this.char_id, str_cell);
            
            % Maybe in a future it will be useful to export and import specific satellites in use
            % str_cell = Ini_Manager.toIniString('index', this.index, str_cell);
            % str_cell = Ini_Manager.toIniString('prn', this.prn, str_cell);
            % str_cell = Ini_Manager.toIniString('system', this.system, str_cell);
        end                
    end
    
    % =========================================================================
    %  LEGACY IMPORT
    % =========================================================================
    methods (Access = 'public')
        function legacyImport(this, state)
            % import from the state variable (saved into the old interface mat file of goGPS)
            % This function import processing settings from another setting object or ini file
            enabled_ss = false(1,this.N_SYS_TOT);
            enabled_ss(state.activeGNSS) = true;
            enabled_ss(6) = enabled_ss(6) || state.use_sbas;
            this.init(enabled_ss);
        end
    end
            
    % =========================================================================
    %  TEST
    % =========================================================================    
    methods (Static, Access = 'public')        
        function test()      
            % test the class
            c = Constellation_Collector(Constellation_Collector.SYS_C);
            c.testInterfaceRoutines();
        end
    end    

end
