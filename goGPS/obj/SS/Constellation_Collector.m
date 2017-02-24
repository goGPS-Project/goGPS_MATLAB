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
        SYS_EXT_NAME = {'GPS', 'GLONASS', 'Galileo', 'BeiDou', 'QZSS', 'SBAS'}; % full name of the constellation
        SYS_NAME     = {'GPS', 'GLO', 'GAL', 'BDS', 'QZS', 'SBS'}; % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        ID_GPS     = 1 % Id of GPS constellation for goGPS internal use
        ID_GLONASS = 2 % Id of GLONASS constellation for goGPS internal use
        ID_GALILEO = 3 % Id of Galileo constellation for goGPS internal use
        ID_BEIDOU  = 4 % Id of BeiDou constellation for goGPS internal use
        ID_QZSS    = 5 % Id of QZSS constellation for goGPS internal use
        ID_SBAS    = 6 % Id of SBAS constellation for goGPS internal use
        SYS_C      = 'GRECJS'; % Array of constellation char ids GPS = 'G', GLONASS = 'R', ...        
    end
    
    properties (SetAccess = private, GetAccess = public)
        list         % struct of objects (to keep the name of the active constellations)
        sys_name     % array of active constellations names (as in the list structure)
        num_id       % array of active constellations numeric id
        sys_c        % array of active constellations char id
        active_list  % logical array of satellite actually active, order: GPS GLO GAL BDS QZS SBS
        n_sys        % number of active constellations
        n_sat        % uint8 array of satellite used per active constellation
        n_sat_tot    % uint8 total teoretical number of satellite available for processing
        
        index        % incremental index of the active satellite system
        prn          % relative id number in the satellite system        
        system       % char id of the constellation per satellite
    end
    
    methods (Access = 'private')        
        function init(this, active_ss_list)
            % Init function for the Constellation_Collector Class
            % SYNTAX: this.init(active_ss_list)
                this.list.GPS = GPS_SS();
                this.list.GLO = GLONASS_SS();
                this.list.GAL = Galileo_SS();
                this.list.BDS = BeiDou_SS();
                this.list.QZS = QZSS_SS();
                this.list.SBS = SBAS_SS();
                this.updateStatus(active_ss_list);
        end
        
        function updateStatus(this, active_ss_list)
            % Update function for the Constellation_Collector Class
            % SYNTAX: this.updateStatus(active_ss_list)            
            this.active_list = active_ss_list;
            if (sum(active_ss_list) > 0)
                this.prn = [];     % relative id number in the satellite system
                this.system = '';  % char id of the constellation per satellite
                
                this.n_sat_tot = 0; % counter for number of satellites                
                this.num_id = [];
                this.sys_c = [];
                this.system = [];
                this.prn = [];
                this.n_sat = [];
                this.n_sat_tot = 0;                
                if active_ss_list(1) % GPS is active
                    this.list.GPS.updateGoIds(this.n_sat_tot);
                    this.list.GPS.enable();
                    this.num_id = [this.num_id this.ID_GPS];
                    this.sys_c = [this.sys_c this.list.GPS.SYS_C];
                    this.system = [this.system char(ones(1, this.list.GPS.N_SAT) * this.list.GPS.SYS_C)];
                    this.prn = [this.prn; this.list.GPS.PRN];
                    this.n_sat = [this.n_sat this.list.GPS.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.list.GPS.N_SAT;
                end
                if active_ss_list(2) % GLONASS is active
                    this.list.GLO.updateGoIds(this.n_sat_tot);
                    this.list.GLO.enable();
                    this.num_id = [this.num_id this.ID_GLONASS];
                    this.sys_c = [this.sys_c this.list.GLO.SYS_C];
                    this.system = [this.system char(ones(1, this.list.GLO.N_SAT) * this.list.GLO.SYS_C)];
                    this.prn = [this.prn; this.list.GLO.PRN];
                    this.n_sat = [this.n_sat this.list.GLO.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.list.GLO.N_SAT;
                end
                if active_ss_list(3) % Galileo is active
                    this.list.GAL.updateGoIds(this.n_sat_tot);
                    this.list.GAL.enable();
                    this.num_id = [this.num_id this.ID_GALILEO];
                    this.sys_c = [this.sys_c this.list.GAL.SYS_C];
                    this.system = [this.system char(ones(1, this.list.GAL.N_SAT) * this.list.GAL.SYS_C)];
                    this.prn = [this.prn; this.list.GAL.PRN];
                    this.n_sat = [this.n_sat this.list.GAL.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.list.GAL.N_SAT;
                end
                if active_ss_list(4) % BeiDou is active
                    this.list.BDS.updateGoIds(this.n_sat_tot);
                    this.list.BDS.enable();
                    this.num_id = [this.num_id this.ID_BEIDOU];
                    this.sys_c = [this.sys_c this.list.BDS.SYS_C];
                    this.system = [this.system char(ones(1, this.list.BDS.N_SAT) * this.list.BDS.SYS_C)];
                    this.prn = [this.prn; this.list.BDS.PRN];
                    this.n_sat = [this.n_sat this.list.BDS.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.list.BDS.N_SAT;
                end
                if active_ss_list(5) % QZSS is active
                    this.list.QZS.updateGoIds(this.n_sat_tot);
                    this.list.QZS.enable();
                    this.num_id = [this.num_id this.ID_QZSS];
                    this.sys_c = [this.sys_c this.list.QZS.SYS_C];
                    this.system = [this.system char(ones(1, this.list.QZS.N_SAT) * this.list.QZS.SYS_C)];
                    this.prn = [this.prn; this.list.QZS.PRN];
                    this.n_sat = [this.n_sat this.list.QZS.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.list.QZS.N_SAT;
                end
                if active_ss_list(6) % SBAS is active (not yet implemented)
                    this.list.SBS.updateGoIds(this.n_sat_tot);
                    this.list.SBS.enable();
                    this.num_id = [this.num_id this.ID_SBAS];
                    this.sys_c = [this.sys_c this.list.SBS.SYS_C];
                    this.system = [this.system char(ones(1, this.list.SBS.N_SAT) * this.list.SBS.SYS_C)];
                    this.prn = [this.prn; this.list.SBS.PRN];
                    this.n_sat = [this.n_sat this.list.SBS.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.list.SBS.N_SAT;
                end
                
                this.index = (1 : this.n_sat_tot)';   % incremental index of the active satellite system
                this.n_sys = numel(this.num_id);
                this.sys_name = this.SYS_NAME(this.num_id);
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
                        active_ss_list = false(1,this.N_SYS_TOT);
                        [~, ids] = intersect('GRECJS',GPS_flag);
                        active_ss_list(ids) = true;
                    else
                        active_ss_list = logical(GPS_flag);
                    end
                case 5,  active_ss_list = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, 0]);
                case 6,  active_ss_list = logical([GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag]);
                otherwise, error(['Initialization of Constellation_Collector failed: ' 10 '   invalid number of parameters in the constructor call']);
            end
            
            % check the size of the array active_list
            if (numel(active_ss_list) < this.N_SYS_TOT)
                tmp = false(this.N_SYS_TOT, 1);
                tmp(1:numel(active_ss_list)) = active_ss_list;
                active_ss_list = tmp;
                clear tmp;
            else
                active_ss_list = active_ss_list(1:this.N_SYS_TOT);
            end
            if active_ss_list(6)
                this.logger.addWarning('usage of SBAS satellite system not yet tested in the latest goGPS version, something may not work correctly!');
            end
            this.init(active_ss_list);
        end
        
        function update(this)
            % update the class, if some constallation have been modified
            % SYNTAX: this.update();
            this.getActive;
        end

        function active_ss_list = getActive(this)
            % get the logical array of satellite actually active, order: GPS GLO GAL BDS QZS SBS -> and update the object if something has changed
            % SYNTAX: active_list = this.getActive();
            active_ss_list = [this.list.GPS.isActive this.list.GLO.isActive this.list.GAL.isActive this.list.BDS.isActive this.list.QZS.isActive this.list.SBS.isActive];
            
            % If some constellation have been activated not in the proper way 
            if (sum(not(this.active_list & active_ss_list)) > 0) 
                this.updateStatus(active_ss_list);
            end
        end

    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================    
    methods (Access = 'public')
        function import(this, settings)
            % This function import processing settings from another setting object or ini file
            if isa(settings, 'Ini_Manager')
                vl = this.logger.getVerbosityLev();
                this.logger.setVerbosityLev(0);
                this.init([0 0 0 0 0 0]);
            	this.logger.setVerbosityLev(vl);
                
                this.list.GPS.import(settings);
                this.list.GLO.import(settings);
                this.list.GAL.import(settings);
                this.list.BDS.import(settings);
                this.list.QZS.import(settings);
                this.list.SBS.import(settings);                
                this.update();
            else
                this.list = repmat(settings.list,1,1);
                this.update();
            end
        end
        
        function str = toString(this, str)
            % Display the satellite system in use            
            if (nargin == 1)
                str = '';
            end
            this.update();
            [~, ids] = intersect(this.sys_c, this.sys_c);
            toString = @(var) regexprep(regexprep(evalc(['disp(var)']), '''   ', ','), '''', '');
            %str = [str ' Constellation in use: ' toString(this.SYS_EXT_NAME(sort(ids)))];
            
            str = this.list.GPS.toString(str);
            str = this.list.GLO.toString(str);
            str = this.list.GAL.toString(str);
            str = this.list.BDS.toString(str);
            str = this.list.QZS.toString(str);
            str = this.list.SBS.toString(str);
            str = [str 10];
        end
        
        function str_cell = export(this, str_cell)
            % Conversion to string of the minimal information needed to reconstruct the this            
            if (nargin == 1)
                str_cell = {};
            end
            this.update();
            str_cell = Ini_Manager.toIniStringComment('Constallations for the processing:', str_cell);
            %str_cell = Ini_Manager.toIniStringComment('insert the constellations as array of char:', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "G" GPS', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "R" GLONASS', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "E" Galileo', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "C" BeiDou', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "J" QZSS', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "S" SBAS (not yet available)', str_cell);
            %str_cell = Ini_Manager.toIniString('active_constellation_ch', this.sys_c, str_cell);
            
            str_cell = this.list.GPS.export(str_cell);
            str_cell = this.list.GLO.export(str_cell);
            str_cell = this.list.GAL.export(str_cell);
            str_cell = this.list.BDS.export(str_cell);
            str_cell = this.list.QZS.export(str_cell);
            str_cell = this.list.SBS.export(str_cell);
            
            % Maybe in a future it will be useful to export and import specific satellites in use
            % str_cell = Ini_Manager.toIniString('index', this.index, str_cell);
            % str_cell = Ini_Manager.toIniString('PRN', this.prn, str_cell);
            % str_cell = Ini_Manager.toIniString('system', this.system, str_cell);
        end                
    end
    
    % =========================================================================
    %  LEGACY IMPORT
    % =========================================================================
    methods (Access = 'public')
        function legacyImport(this, state)
            % Import from the state variable (saved into the old interface mat file of goGPS)
            % This function import processing settings from another setting object or ini file
            % SYNTAX: obj.legacyImport(state)
            active_ss_list = false(1,this.N_SYS_TOT);
            active_ss_list(state.activeGNSS) = true;
            active_ss_list(6) = active_ss_list(6) || state.use_sbas;
            this.init(active_ss_list);
            if isfield(state,'activeFreq')
                this.list.GPS.setActiveFrequencies(state.activeFreq);
                this.list.GLO.setActiveFrequencies(state.activeFreq);
                this.list.GAL.setActiveFrequencies(state.activeFreq);
                this.list.BDS.setActiveFrequencies(state.activeFreq);
                this.list.QZS.setActiveFrequencies(state.activeFreq);
            end
        end
    end
            
    % =========================================================================
    %  TEST
    % =========================================================================    
    methods (Static, Access = 'public')        
        function test()
            % Test the class
            % SYNTAX: Constellation_Collector.test()
            c = Constellation_Collector(Constellation_Collector.SYS_C);
            c.testInterfaceRoutines();
            c = Constellation_Collector([0 1 1 0 0 0]);
            c.testInterfaceRoutines();            
        end
    end    

end
