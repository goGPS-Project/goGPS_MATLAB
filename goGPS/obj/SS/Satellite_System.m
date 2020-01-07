%   CLASS Satellite_System
% =========================================================================
%
% DESCRIPTION
%   Abstract class with basic properties and methods that a satellite
%   system must have
%

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b5 Happy 2020
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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

classdef Satellite_System < Settings_Interface
    properties (Constant, Abstract)
        SYS_EXT_NAME  % full name of the constellation
        SYS_NAME      % 3 characters name of the constellation, this "short name" is used as fields of the property list (struct) to identify a constellation
        SYS_C         % Satellite system (ss) character id

        F             % System frequencies as struct
        F_VEC         % Array of supported frequencies
        L_VEC         % Array of the corresponding wavelength (pre-computed for performance reason)

        N_SAT         % Maximum number of satellite in the constellation
        PRN           % Satellites id numbers as defined in the constellation

        % Structure of orbital parameters (ellipsoid, GM, OMEGA_EARTH_DOT)
        ORBITAL_P
        %   .GM
        %   .OMEGAE_DOT
        %   .J2 (only for GLONASS)
        %   .ELL.a      semi-major axis
        %   .ELL.F      flattening
        %   .ELL.e      eccentricity
        %   .ELL.e2     eccentricity^2
        
        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB;  % last letter of the observation code e.g. IRNSS: C5A - C5B - C5C - C5X
        CODE_RIN3_2BAND;  % id for the freq as stored in F_VEC e.g. IRNSS: L5 -> C5A, S -> C9A
    end

    properties (GetAccess = 'public', SetAccess = 'protected')
        % Flag array containing the list of active frequencies
        flag_f = [];

        % flag defining the status of activation of the constellation
        flag_enable = false;

        % Satellites unique id numbers in goGPS
        go_ids
    end

    methods (Access = public)
        function iono_free = getIonoFree(this, f_id)
            % Init iono parameters: iono-free combination is computed with the first two carriers in F_VEC (use f_id to change the frequencies to use)
            % 
            % SYNTAX
            %   iono_free = ss_obj.getIonoFree(<f_id>)
            %
            % Structure of the iono free combination parameters (alpha 1, alpha 2, T ed N)
            % iono_free
            %   .alpha1
            %   .alpha2
            %   .T
            %   .N
            if nargin < 2
                f_id = [1 2];
            end
            iono_free.alpha1 = this.F_VEC(f_id(1)) .^ 2 ./ (this.F_VEC(f_id(1)) .^ 2 - this.F_VEC(f_id(2)) .^ 2);
            iono_free.alpha2 = this.F_VEC(f_id(2)) .^ 2 ./ (this.F_VEC(f_id(1)) .^ 2 - this.F_VEC(f_id(2)) .^ 2);
            gcd_f = gcd(this.F_VEC(f_id(1)),this.F_VEC(f_id(2)));
            iono_free.T = this.F_VEC(f_id(1))/gcd_f;
            iono_free.N = this.F_VEC(f_id(2))/gcd_f;
        end
    end

    methods
        function this = Satellite_System(offset)
            % Creator
            %
            % SYNTAX
            %   Satellite_System(offset)
            if (nargin == 0)
                offset = 0;
            end
            this.flag_f = false(1, size(this.F_VEC, 2));
            this.updateGoIds(offset);
            this.setActiveFrequencies([1 0 0]);
        end

        function setActiveFrequencies(this, flag_freq)
            % Set the frequencies that are going to be use
            %
            % SYNTAX
            %   setActiveFrequencies(flag_freq)
            len = min(length(flag_freq), length(this.flag_f));
            flags =  false(numel(this.flag_f),1);
            flags(1:len) = flag_freq(1:len);
            this.flag_f = flags;
        end

        function updateGoIds(this, offset)
            % Update the satellites unique id numbers in goGPS
            %
            % SYNTAX
            %   updateGoIds(<offset>)
            if nargin == 1
                offset = 0;
            end
            this.go_ids = offset + (1 : this.N_SAT);
        end

        function code = getPrCodes(this, freq_num)
            % get the list of codes containing pseudo range data
            %
            % SYNTAX
            %   code = getPrCodes(this, freq_num)
            code_rin3_avail = this.CODE_RIN3_ATTRIB{freq_num};
            code_rin3_avail(code_rin3_avail == 'N') = []; % remove codeless observations
            code = char([ones(numel(code_rin3_avail), 1) * 'C' ones(numel(code_rin3_avail), 1) * this.CODE_RIN3_2BAND(freq_num) code_rin3_avail']);
        end
        
        function code = getPhCodes(this, freq_num)
            % get the list of codes containing phase data
            %
            % SYNTAX
            %   code = getPhCodes(this, freq_num)
            code_rin3_avail = this.CODE_RIN3_ATTRIB{freq_num};
            code = char([ones(numel(code_rin3_avail), 1) * 'L' ones(numel(code_rin3_avail), 1) * this.CODE_RIN3_2BAND(freq_num) code_rin3_avail']);
        end
        
        function code = getDopCodes(this, freq_num)
            % get the list of codes containing doppler
            %
            % SYNTAX
            %   code = getDopCodes(this, freq_num)
            code_rin3_avail = this.CODE_RIN3_ATTRIB{freq_num};
            code = char([ones(numel(code_rin3_avail), 1) * 'D' ones(numel(code_rin3_avail), 1) * this.CODE_RIN3_2BAND(freq_num) code_rin3_avail']);
        end
        
        function code = getSnrCodes(this, freq_num)
            % get the list of codes containing snr data
            %
            % SYNTAX
            %   code = getSnrCodes(this, freq_num)
            code_rin3_avail = this.CODE_RIN3_ATTRIB{freq_num};
            code = char([ones(numel(code_rin3_avail), 1) * 'S' ones(numel(code_rin3_avail), 1) * this.CODE_RIN3_2BAND(freq_num) code_rin3_avail']);
        end
        
        function offset = getOffset(this)
            % Get offset of the go_ids
            %
            % SYNTAX
            %   offset = getOffset(this)
            if ~isempty(this.go_ids)
                offset = this.go_ids(1) - 1;
            else
                offset = 0;
            end
        end

        function id = getFirstId(this)
            % get the first goGPS id -> if constellation is inactive
            %
            % SYNTAX
            %   id = getFirstId(this)
            if this.flag_enable % this.isActive()
                id = this.go_ids(1);
            else
                id = 0;
            end
        end

        function name = getFreqName(this)
            % Get the name of the frequencies for the current constellation
            %
            % SYNTAX
            %   name = this.getFreqNames();
            name = fieldnames(this.F);
        end
        
        function setFlagF(this, idx, values)
            this.flag_f(idx) = values;
        end

        function flag = isActive(this)
            % Get the status of activation of the constellation
            %
            % SYNTAX
            %   flag = isActive(this)
            flag = this.flag_enable;
        end

        function enable(this, status)
            % Enable this satellite system
            % Warning: when the obj is used within Constellation Collectore  use activateXXXX deactivateXXXX functions
            %
            % SYNTAX
            %   this.enable(<status>);
            if (nargin == 2)
                this.flag_enable = logical(status);
                if isempty(this.flag_enable)
                    this.flag_enable = false;
                end
            else
                this.flag_enable = true;
            end
        end

        function disable(this)
            % Disable this satellite system
            % Warning: when the obj is used within Constellation Collectore
            % use activateXXXX deactivateXXXX functions
            %
            % SYNTAX
            %   this.disable();
            this.flag_enable = false;
        end
    end

    methods (Abstract)
        copy = getCopy(this);
    end

    % =========================================================================
    %  INTERFACE REQUIREMENTS
    % =========================================================================
    methods
        function import(this, settings)
            % This function import Satellite System (only) settings from another setting object
            %
            % SYNTAX
            %   this.import(settings)
            if isa(settings, 'Ini_Manager')
                this.enable(settings.getData(sprintf('%s_is_active', this.SYS_NAME)));
                name = this.getFreqName();
                tmp = true(numel(name),1);
                for i = 1 : numel(name)
                    tmp_str = logical(settings.getData(sprintf('%s_%s', this.SYS_NAME, char(name{i}))));
                    if ~isempty(tmp_str)
                        tmp(i) = tmp_str;
                    end
                end
                this.flag_f = tmp;
            else
                this.flag_f = settings.flag_f;
                this.enable(settings.flag_enable);
                this.go_ids = settings.go_ids;
            end
        end

        function str = toString(this, str)
            % Display the satellite system in use
            %
            % SYNTAX
            %   str = toString(this, str)
            if (nargin == 1)
                str = '';
            end
            if (this.isActive())
                name = this.getFreqName();
                str = [str sprintf(' Satellite system "%s" active\n  - using frequencies: %s\n', this.SYS_EXT_NAME, sprintf('"%s" ',  strCell2Str(name(this.flag_f),', ')))];
            else
                str = [str sprintf(' Satellite system "%s" is inactive\n', this.SYS_EXT_NAME)];
            end
        end

        function str_cell = export(this, str_cell)
            % Conversion to string ini format of the minimal information needed to reconstruct the this
            %
            % SYNTAX
            %   str = toString(this, str)
            if (nargin == 1)
                str_cell = {};
            end
            name = this.getFreqName();
            str_cell = Ini_Manager.toIniStringComment(sprintf('%s satellite system', this.SYS_EXT_NAME), str_cell);
            str_cell = Ini_Manager.toIniString(sprintf('%s_is_active', this.SYS_NAME), this.isActive(), str_cell);
            str_cell = Ini_Manager.toIniStringComment('Frequencies to be used when this constellation is active', str_cell);
            for f = 1 : numel(name)
                str_cell = Ini_Manager.toIniString(sprintf('%s_%s', this.SYS_NAME, name{f}), this.isActive() && this.flag_f(f), str_cell);
            end
            str_cell = Ini_Manager.toIniStringNewLine(str_cell);
        end
    end

    % =========================================================================
    %  LEGACY IMPORT
    % =========================================================================
    methods (Access = 'public')
        function legacyImport(this, state)
            % import from the state variable (saved into the old interface mat file of goGPS)
            %
            % SYNTAX
            %   this.legacyImport(state)
            try
                this.p_mode = this.gui2mode(state.mode, state.nav_mon, state.kalman_ls, state.code_dd_sa);
            catch ex
                this.log.addWarning(['Legacy import "Processing mode" failed - ', ex.message])
            end
        end
    end
end
