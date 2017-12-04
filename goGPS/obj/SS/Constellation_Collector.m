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
%     __ _ ___ / __| _ | __
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 0.5.2 beta 1
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
    properties (Constant, GetAccess = public)
        N_SYS_TOT = 7; % Max number of available satellite systems
        SYS_EXT_NAME = {'GPS', 'GLONASS', 'Galileo', 'QZSS', 'BeiDou', 'IRNSS', 'SBAS'}; % full name of the constellation
        SYS_NAME     = {'GPS', 'GLO', 'GAL', 'QZS', 'BDS', 'IRN', 'SBS'}; % 3 characters name of the constellation
        ID_GPS       = 1 % Id of GPS constellation for goGPS internal use
        ID_GLONASS   = 2 % Id of GLONASS constellation for goGPS internal use
        ID_GALILEO   = 3 % Id of Galileo constellation for goGPS internal use
        ID_QZSS      = 4 % Id of QZSS constellation for goGPS internal use
        ID_BEIDOU    = 5 % Id of BeiDou constellation for goGPS internal use
        ID_IRNSS     = 6 % Id of IRNSS constellation for goGPS internal use
        ID_SBAS      = 7 % Id of SBAS constellation for goGPS internal use
        SYS_C        = 'GREJCIS'; % Array of constellation char ids GPS = 'G', GLONASS = 'R', ...
    end

    properties (SetAccess = public, GetAccess = public)
        gps = GPS_SS();      % GPS parameters
        glo = GLONASS_SS();  % GLONASS parameters
        gal = Galileo_SS();  % Galileo parameters
        bds = BeiDou_SS();   % BeiDou parameters
        qzs = QZSS_SS();     % QZSS parameters
        irn = IRNSS_SS();    % IRNSS parameters
        sbs = SBAS_SS();     % sbs parameters

        sys_name     % array of active constellations names (as in the list structure)
        num_id       % array of active constellations numeric id
        sys_c        % array of active constellations char id
        active_list  % logical array of satellite actually active, order: GPS GLO GAL BDS QZS SBS
        n_sys        % number of active constellations
        n_sat        % uint8 array of satellite used per active constellation
        n_sat_tot    % uint8 total teoretical number of satellite available for processing
        n_obs_tot    % uint16 total number of observables

        index        % incremental index of the active satellite system
        prn          % relative id number in the satellite system
        system       % char id of the constellation per satellite
    end

    methods (Access = 'private')
        function init(this, active_list)
            % Init function for the Constellation_Collector Class
            % SYNTAX: this.init(active_list)

            this.updateStatus(active_list);
        end

        function updateStatus(this, active_list)
            % Update function for the Constellation_Collector Class
            % SYNTAX: this.updateStatus(active_list)
            this.active_list = active_list;
            if (sum(active_list) > 0)
                this.prn = [];     % relative id number in the satellite system
                this.system = '';  % char id of the constellation per satellite

                this.n_sat_tot = 0; % counter for number of satellites
                this.num_id = [];
                this.sys_c = [];
                this.system = [];
                this.prn = [];
                this.n_sat = [];
                this.n_sat_tot = 0;
                this.n_obs_tot = 0;
                if active_list(1) % GPS is active
                    this.gps.updateGoIds(this.n_sat_tot);
                    this.gps.enable();
                    this.num_id = [this.num_id this.ID_GPS];
                    this.sys_c = [this.sys_c this.gps.SYS_C];
                    this.system = [this.system char(ones(1, this.gps.N_SAT) * this.gps.SYS_C)];
                    this.prn = [this.prn; this.gps.PRN];
                    this.n_sat = [this.n_sat this.gps.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.gps.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.gps.N_SAT * sum(this.gps.flag_f);
                else
                    this.gps.disable();
                end
                if active_list(2) % GLONASS is active
                    this.glo.updateGoIds(this.n_sat_tot);
                    this.glo.enable();
                    this.num_id = [this.num_id this.ID_GLONASS];
                    this.sys_c = [this.sys_c this.glo.SYS_C];
                    this.system = [this.system char(ones(1, this.glo.N_SAT) * this.glo.SYS_C)];
                    this.prn = [this.prn; this.glo.PRN];
                    this.n_sat = [this.n_sat this.glo.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.glo.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.glo.N_SAT * sum(this.glo.flag_f);
                else
                    this.glo.disable();
                end
                if active_list(3) % Galileo is active
                    this.gal.updateGoIds(this.n_sat_tot);
                    this.gal.enable();
                    this.num_id = [this.num_id this.ID_GALILEO];
                    this.sys_c = [this.sys_c this.gal.SYS_C];
                    this.system = [this.system char(ones(1, this.gal.N_SAT) * this.gal.SYS_C)];
                    this.prn = [this.prn; this.gal.PRN];
                    this.n_sat = [this.n_sat this.gal.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.gal.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.gal.N_SAT * sum(this.gal.flag_f);
                else
                    this.gal.disable();
                end
                if active_list(4) % QZSS is active
                    this.qzs.updateGoIds(this.n_sat_tot);
                    this.qzs.enable();
                    this.num_id = [this.num_id this.ID_QZSS];
                    this.sys_c = [this.sys_c this.qzs.SYS_C];
                    this.system = [this.system char(ones(1, this.qzs.N_SAT) * this.qzs.SYS_C)];
                    this.prn = [this.prn; this.qzs.PRN];
                    this.n_sat = [this.n_sat this.qzs.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.qzs.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.qzs.N_SAT * sum(this.qzs.flag_f);
                else
                    this.qzs.disable();
                end
                if active_list(5) % BeiDou is active
                    this.bds.updateGoIds(this.n_sat_tot);
                    this.bds.enable();
                    this.num_id = [this.num_id this.ID_BEIDOU];
                    this.sys_c = [this.sys_c this.bds.SYS_C];
                    this.system = [this.system char(ones(1, this.bds.N_SAT) * this.bds.SYS_C)];
                    this.prn = [this.prn; this.bds.PRN];
                    this.n_sat = [this.n_sat this.bds.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.bds.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.bds.N_SAT * sum(this.bds.flag_f);
                else
                    this.bds.disable();
                end
                if active_list(6) % IRNSS is active
                    this.irn.updateGoIds(this.n_sat_tot);
                    this.irn.enable();
                    this.num_id = [this.num_id this.ID_IRNSS];
                    this.sys_c = [this.sys_c this.irn.SYS_C];
                    this.system = [this.system char(ones(1, this.irn.N_SAT) * this.irn.SYS_C)];
                    this.prn = [this.prn; this.irn.PRN];
                    this.n_sat = [this.n_sat this.irn.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.irn.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.irn.N_SAT * sum(this.irn.flag_f);
                else
                    this.irn.disable();
                end
                
                if active_list(7) % SBAS is active (not yet implemented)
                    this.sbs.updateGoIds(this.n_sat_tot);
                    this.sbs.enable();
                    this.num_id = [this.num_id this.ID_SBAS];
                    this.sys_c = [this.sys_c this.sbs.SYS_C];
                    this.system = [this.system char(ones(1, this.sbs.N_SAT) * this.sbs.SYS_C)];
                    this.prn = [this.prn; this.sbs.PRN];
                    this.n_sat = [this.n_sat this.sbs.N_SAT];
                    this.n_sat_tot = this.n_sat_tot + this.sbs.N_SAT;
                    this.n_obs_tot = this.n_obs_tot + this.sbs.N_SAT * sum(this.sbs.flag_f);
                else
                    this.sbs.disable();
                end

                this.index = (1 : this.n_sat_tot)';   % incremental index of the active satellite system
                this.n_sys = numel(this.num_id);
                this.sys_name = this.SYS_NAME(this.num_id);
            else
                this.log.addError('No satellite system selected -> Enabling GPS');
                ss = false(this.N_SYS_TOT,1); ss(1) = true;
                this.init(ss);
            end
        end
    end

    methods
        function this = Constellation_Collector(GPS_flag, GLO_flag, GAL_flag, QZS_flag, BDS_flag, IRN_flag, SBS_flag)
            % Constructor - parameters: [GPS_flag, GLO_flag, GAL_flag, QZS_flag, BDS_flag, IRN_flag, SBS_flag]
            % SYNTAX:
            %   cc = Constellation_Collector('GRECJS'); % use the array of constellation char ids
            %   cc = Constellation_Collector(GPS_flag, GLO_flag, GAL_flag, QZS_flag, BDS_flag, IRN_flag, SBS_flag);
            %   cc = Constellation_Collector([GPS_flag, GLO_flag, GAL_flag, QZS_flag, BDS_flag, IRN_flag, SBS_flag]);
            %
            % INPUT:
            %   single logical array whose elements are:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   IRN_flag = boolean flag for enabling/disabling IRNSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   object handler
            %
            % DESCRIPTION:
            %   Multi-constellation set-up.

            this.initLogger();

            % this f***ing software called MATLAB does not create new
            % constellation objects unless I do this
            this.gps = this.gps.getCopy();
            this.glo = this.glo.getCopy();
            this.gal = this.gal.getCopy();
            this.qzs = this.qzs.getCopy();
            this.bds = this.bds.getCopy();
            this.irn = this.irn.getCopy();
            this.sbs = this.sbs.getCopy();

            % Manually manage overloading
            switch nargin
                case 0
                    active_list = false(1,this.N_SYS_TOT);
                case 1
                    if (ischar(GPS_flag))
                        active_list = false(1,this.N_SYS_TOT);
                        [~, ids] = intersect('GREJCIS',GPS_flag);
                        active_list(ids) = true;
                    else
                        active_list = logical(GPS_flag);
                    end
                case 6,  active_list = logical([GPS_flag, GLO_flag, GAL_flag, QZS_flag, BDS_flag, IRN_flag, 0]);
                case 7,  active_list = logical([GPS_flag, GLO_flag, GAL_flag, QZS_flag, BDS_flag, IRN_flag, SBS_flag]);
                otherwise, error(['Initialization of Constellation_Collector failed: ' 10 '   invalid number of parameters in the constructor call']);
            end

            % check the size of the array active_list
            if (numel(active_list) < this.N_SYS_TOT)
                tmp = false(this.N_SYS_TOT, 1);
                tmp(1:numel(active_list)) = active_list;
                active_list = tmp;
                clear tmp;
            else
                active_list = active_list(1:this.N_SYS_TOT);
            end
            if active_list(7)
                this.log.addWarning('usage of SBAS satellite system not yet tested in the latest goGPS version, something may not work correctly!');
            end
            this.init(active_list);
        end

        function update(this)
            % update the class, if some constallation have been modified
            % SYNTAX: this.update();
            this.getActive();
        end

        function active_list = getActive(this)
            % get the logical array of satellite actually active, order: GPS GLO GAL QZS BDS SBS -> and update the object if something has changed
            % SYNTAX: active_list = this.getActive();
            active_list = [this.gps.isActive this.glo.isActive this.gal.isActive this.qzs.isActive this.bds.isActive this.irn.isActive this.sbs.isActive];

            % If some constellation have been activated not in the proper way
            if (sum(not(this.active_list(:) & active_list(:))) > 0)
                this.updateStatus(active_list);
            end
        end

        function ant_mod = getAntennaId(this, goid)
            % Equivalent to the old function sat_antenna_ID
            % return the "name" of the antennas of each active satellite
            % SYNTAX:
            %   ant_mod = getAntennaId(this)
            % SEE ALSO:
            %   sat_antenna_ID
            prn = this.prn;
            prn(this.system == 'J') = prn(this.system == 'J') - 192;
            ant_mod = mat2cell([this.system' reshape(sprintf('%02d', prn), 2, this.n_sat_tot)'], ones(this.n_sat_tot,1));
            if nargin > 1
                ant_mod = ant_mod{goid},
            end
        end
        
        function lambda = getLegacyLambda(this, eph, sp3)
            % Get lambda matrix [n_sat x 7], used in a lot of old goGPS functions
            lambda = zeros(this.n_sat_tot, 7);
            % Get satellites with orbits
            if (isempty(sp3))
                sat_ok = unique(eph(30,:))';
            else
                sat_ok = sp3.avail(sp3.avail > 0);
            end
            if this.glo.isActive()
                % Transmission channel for glonass given the prn
                channel_id = zeros(size(this.glo.go_ids));
                for i = 1 : numel(this.glo.go_ids)
                    pos = find(eph(30,:) == this.glo.go_ids(i), 1, 'first');
                    channel_id(i) = eph(15, pos) + 8; % Add 8 cause the GLONASS channel start from -7
                end
            end
            
            for s = 1 : numel(this.sys_name)
                pos = this.(char((this.sys_name{s} + 32))).go_ids();
                iono_free = this.(char((this.sys_name{s} + 32))).getIonoFree();
                if this.sys_c(s) == 'R' % if glonass
                    lambda(pos, 1) = this.(char((this.sys_name{s} + 32))).L_VEC(channel_id,1);
                    lambda(pos, 2) = this.(char((this.sys_name{s} + 32))).L_VEC(channel_id,2);
                else
                    lambda(pos, 1) = this.(char((this.sys_name{s} + 32))).L_VEC(:,1);
                    lambda(pos, 2) = this.(char((this.sys_name{s} + 32))).L_VEC(:,2);
                end
                lambda(pos, 3) = (iono_free.T / (iono_free.T^2 - iono_free.N^2)) * lambda(pos,1);
                lambda(pos, 4) = iono_free.alpha1;
                lambda(pos, 5) = iono_free.alpha2;
                lambda(pos, 6) = iono_free.T;
                lambda(pos, 7) = iono_free.N;
            end
            sat_ko = true(this.n_sat_tot, 1); sat_ko(sat_ok) = false;
            lambda(sat_ko, :) = 0;
        end
        function is_ref = isRefFrequency(this,sys,freq)
            %DESCRIPTION: tell if the frequency is used as reference
            %frequency (also in combination with other one) for the
            %computation of sta clocks
            is_ref = false;
            switch(sys)
                case 'G'
                    if freq == 1 | freq == 2
                        is_ref = true;
                    end
                case 'R'
                    if freq == 1 | freq == 2
                        is_ref = true;
                    end
                case 'E'
                    if freq == 1 | freq == 5
                        is_ref = true;
                    end
                case 'J'
                    if freq == 1 | freq == 2
                        is_ref = true;
                    end
                case 'C'
                    if freq == 2 | freq == 6
                        is_ref = true;
                    end
                case 'I'
                    if freq == 5 | freq == 9
                        is_ref = true;
                    end
                case 'S'
                    if freq == 1 | freq == 5
                        is_ref = true;
                    end
            end
        end
        function ref_dcb = getRefDCB(this, goid)
            % consider to move change them on the base of the source of
            % precise orbit
            sys = this.system(goid);
            switch sys
                case 'G'
                    ref_dcb = ['C1WC2W'; 'C1PC2P']; %b
                case 'R'
                    ref_dcb = 'C1PC2P';
                case 'E'
                    ref_dcb = ['C1CC5Q'; 'C1XC5X'];
                case 'C'
                    ref_dcb = 'C2IC7I';
                case 'J'
                    ref_dcb = 'C1XC2X'; % to be tested
                case 'I'
                    ref_dcb = 'C5AC9A';% to be tested
                otherwise
            end
        end

    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods (Access = 'public')
        function import(this, state)
            % This function import processing settings from another setting object or ini file
            if isa(state, 'Ini_Manager')
                vl = this.log.getVerbosityLev();
                this.log.setVerbosityLev(0);
                this.init([0 0 0 0 0 0 0]);
            	this.log.setVerbosityLev(vl);

                this.gps.import(state);
                this.glo.import(state);
                this.gal.import(state);
                this.qzs.import(state);
                this.bds.import(state);
                this.irn.import(state);
                this.sbs.import(state);
                this.update();
            else
                this.gps = state.gps.getCopy();
                this.glo = state.glo.getCopy();
                this.gal = state.gal.getCopy();
                this.qzs = state.qzs.getCopy();
                this.bds = state.bds.getCopy();
                this.irn = state.qzs.getCopy();
                this.sbs = state.sbs.getCopy();
                this.update();
            end
        end

        function str = toString(this, str)
            % Display the satellite system in use
            if (nargin == 1)
                str = '';
            end
            this.update();
            %[~, ids] = intersect(this.sys_c, this.sys_c);
            %toString = @(var) regexprep(regexprep(evalc(['disp(var)']), '''   ', ','), '''', '');
            %str = [str ' Constellation in use: ' toString(this.SYS_EXT_NAME(sort(ids)))];

            str = this.gps.toString(str);
            str = this.glo.toString(str);
            str = this.gal.toString(str);
            str = this.qzs.toString(str);
            str = this.bds.toString(str);
            str = this.irn.toString(str);
            str = this.sbs.toString(str);
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
            %str_cell = Ini_Manager.toIniStringComment(' - "J" QZSS', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "C" BeiDou', str_cell);
            %str_cell = Ini_Manager.toIniStringComment(' - "S" SBAS (not yet available)', str_cell);
            %str_cell = Ini_Manager.toIniString('active_constellation_ch', this.sys_c, str_cell);

            str_cell = this.gps.export(str_cell);
            str_cell = this.glo.export(str_cell);
            str_cell = this.gal.export(str_cell);
            str_cell = this.qzs.export(str_cell);
            str_cell = this.bds.export(str_cell);
            str_cell = this.irn.export(str_cell);
            str_cell = this.sbs.export(str_cell);

            % Maybe in a future it will be useful to export and import specific satellites in use
            % str_cell = Ini_Manager.toIniString('index', this.index, str_cell);
            % str_cell = Ini_Manager.toIniString('PRN', this.prn, str_cell);
            % str_cell = Ini_Manager.toIniString('system', this.system, str_cell);
        end
    end

    % =========================================================================
    %  STATUS MODIFIER
    % =========================================================================
    methods (Access = 'public')
        function setActive(this, active_list)
            this.gps.enable(active_list(1));
            this.glo.enable(active_list(2));
            this.gal.enable(active_list(3));
            this.qzs.enable(active_list(4));
            this.bds.enable(active_list(5));
            this.irn.enable(active_list(6));
            this.sbs.enable(active_list(7));
        end
        
        function activateGPS(this)
            % Activate the use of GPS
            this.gps.enable();
            this.update();
        end

        function deactivateGPS(this)
            % De-activate the use of GPS
            this.gps.disable();
            this.update();
        end

        function activateGLONASS(this)
            % Activate the use of GLONASS
            this.glo.enable();
            this.update();
        end

        function deactivateGLONASS(this)
            % De-activate the use of GLONASS
            this.glo.disable();
            this.update();
        end

        function activateGalileo(this)
            % Activate the use of Galileo
            this.gal.enable();
            this.update();
        end

        function deactivateGalileo(this)
            % De-activate the use of Galileo
            this.gal.disable();
            this.update();
        end

        function activateBeiDou(this)
            % Activate the use of BeiDou
            this.bds.enable();
            this.update();
        end

        function deactivateBeiDou(this)
            % De-activate the use of BeiDou
            this.bds.disable();
            this.update();
        end

        function activateQZSS(this)
            % Activate the use of QZSS
            this.qzs.enable();
            this.update();
        end

        function deactivateQZSS(this)
            % De-activate the use of QZSS
            this.qzs.disable();
            this.update();
        end

        function activateIRNSS(this)
            % Activate the use of IRNSS
            this.qzs.enable();
            this.update();
        end

        function deactivateIRNSS(this)
            % De-activate the use of IRNSS
            this.irn.disable();
            this.update();
        end
        
        function activateSBAS(this)
            % Activate the use of SBAS
            this.gps.enable();
            this.update();
        end

        function deactivateSBAS(this)
            % De-activate the use of SBAS
            this.sbs.disable();
            this.update();
        end
    end

    % =========================================================================
    %  STATUS CHECKER
    % =========================================================================
    methods (Access = 'public')
        function isActive = isGpsActive(this)
            % get the activation status for GPS
            isActive = this.gps.isActive();
        end

        function isActive = isGloActive(this)
            % get the activation status for GLONASS
            isActive = this.glo.isActive();
        end

        function isActive = isGalActive(this)
            % get the activation status for Galileo
            isActive = this.gal.isActive();
        end

        function isActive = isBdsActive(this)
            % get the activation status for Beidou
            isActive = this.bds.isActive();
        end

        function isActive = isQzsActive(this)
            % get the activation status for QZSS
            isActive = this.qzs.isActive();
        end

        function isActive = isIrnActive(this)
            % get the activation status for IRNSS
            isActive = this.irn.isActive();
        end

        function isActive = isSbsActive(this)
            % get the activation status for SBAS
            isActive = this.sbs.isActive();
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
            active_lst = false(1,this.N_SYS_TOT);
            active_lst(state.activeGNSS) = true;
            active_lst(6) = active_lst(6) || state.use_sbas;
            this.init(active_lst);
            if isfield(state,'activeFreq')
                this.gps.setActiveFrequencies(state.activeFreq);
                this.glo.setActiveFrequencies(state.activeFreq);
                this.gal.setActiveFrequencies(state.activeFreq);
                this.qzs.setActiveFrequencies(state.activeFreq);
                this.bds.setActiveFrequencies(state.activeFreq);
                this.irn.setActiveFrequencies(state.activeFreq);
                this.sbs.setActiveFrequencies(state.activeFreq);
            end
        end
    end

    % =========================================================================
    %  GETTERS
    % =========================================================================
    methods (Access = 'public')
        function n_sat = getNumSat(this)
            % return the total number of enabled satellites
            % SYNTAX: n_sat = this.getNumSat();
            n_sat = this.n_sat_tot;
        end
        
        function n_sat = getMaxObsSat(this)
            % return the total maximum number of observables
            % SYNTAX: n_sat = this.getObsSat();
            n_sat = this.n_obs_tot;
        end

        function n_sys = getNumSys(this)
            % return the total number of enabled satellites
            % SYNTAX: n_sat = this.getNumSat();
            n_sys = this.n_sys;
        end

        function gps = getGPS(this)
            % return the GPS satellite system object
            gps = handle(this.gps);
        end

        function gps = getGLONASS(this)
            % return the GLONASS satellite system object
            gps = handle(this.glo);
        end

        function gps = getGalileo(this)
            % return the Galileo satellite system object
            gps = handle(this.gal);
        end

        function gps = getBeiDou(this)
            % return the BeiDou satellite system object
            gps = handle(this.bds);
        end

        function gps = getQZSS(this)
            % return the QZSS satellite system object
            gps = handle(this.qzs);
        end

        function gps = getIRNSS(this)
            % return the IRNSS satellite system object
            gps = handle(this.irn);
        end
        
        function gps = getSBAS(this)
            % return the SBAS satellite system object
            gps = handle(this.sbs);
        end
        function sys = getSys(this, sys)
            % return the system based on is identified char
            switch sys
                    case 'G'
                        sys = this.getGPS();
                    case 'R'
                        sys = this.getGLONASS();
                    case 'E'
                        sys = this.getGalileo();
                    case 'C'
                        sys = this.getBeiDou();
                    case 'J'
                        sys = this.getQZSS();
                    case 'I'
                        sys = this.getIRNSS();
                    case 'S'
                        sys = this.getSBAS();
            end
        end
        function index = getIndex(this,sys,prn)
            % get progressive index for giveng system and prns
            index = zeros(size(prn));
            for i = 1:length(prn)
                p = prn(i);
                index(i) = find(this.system == sys & [this.prn == p]');
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
            c = Constellation_Collector([0 1 1 0 0 1 0]);
            c.testInterfaceRoutines();
        end
    end

end
