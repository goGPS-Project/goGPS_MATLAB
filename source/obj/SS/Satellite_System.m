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
%    |___/                    v 1.0
%
%--------------------------------------------------------------------------
%  Copyright (C) 2023 Geomatics Research & Development srl (GReD)
%  Written by:        Andrea Gatti
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

        ORBITAL_INC    % Orbital inclination
        ORBITAL_RADIUS % Orbital radius
        
        % CODE2DATA ftp://igs.org/pub/data/format/rinex303.pdf
        CODE_RIN3_ATTRIB; % last letter of the observation code e.g. IRNSS: C5A - C5B - C5C - C5X
        CODE_RIN3_2BAND;  % id for the freq as stored in F_VEC e.g. IRNSS: L5 -> C5A, S -> C9A        
    end

    properties (GetAccess = 'public', SetAccess = 'protected')
        % Flag array containing the list of active frequencies
        flag_f = [];

        % flag defining the status of activation of the constellation
        flag_enable = false;
        weight = 1;

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
        
        function [mask_north, mask_south, out_mask] = getPolarMask(this, rec_lat, rec_lon, offset, step)
            % Get the mask of the two palar caps due to the orbital inclination of the satellites
            %
            % INPUT
            %   rec_lat   receiver latitude  [rad]
            %   rec_lon   receiver_longitude [rad]
            %   <offset>  add this to the orbit inclination [deg]
            %   <step>    [el az] grid step [deg] default = 0.5
            %
            % OUTPUT
            %   mask_north     [n x 2] az, el of the sky limits (north)
            %   mask_south     [n x 2] az, el of the sky limits (south)
            %
            % SYNTAX
            %    [mask_north, mask_south, out_mask] = getPolarMask(this, rec_lat, rec_lon, <offset = 0>, <step>)
            if nargin < 4 || isempty(offset)
                offset = 0;
            end
            
            % Elevation resolution
            if nargin < 5 && nargout == 3
                az_step = 0.5;
                el_step = 0.5;
            else
                el_step = step(1);
                az_step=  step(end);
            end
            if nargout == 3
                az_grid = linspace(0, 360-az_step, round(360/az_step))' + az_step/2;
                out_mask = false(round(90/el_step), numel(az_grid));
            end
            
            [mask_north, mask_south] = Satellite_System.generatePolarMask((this.ORBITAL_INC + offset) / 180 * pi, this.ORBITAL_RADIUS, rec_lat, rec_lon);        
            if nargout == 3
                mask_tmp = mask_north;
                for m = 1:2
                    if not(isempty(mask_tmp))
                        mask_tmp(:,1) = mod(mask_tmp(:,1)+pi,2*pi)-pi; % mask az north from -pi to pi
                        id_pos = find(mask_tmp(:,1) >= 0);
                        mask_tmp(:,1) = mod(mask_tmp(:,1),2*pi); % mask az north from 0 to 2*pi
                        
                        [~, id_min] = min(abs(mask_tmp(id_pos,1)));
                        
                        id_min = id_pos(id_min); % norhtern point in a polar plot
                        % making az monothonic
                        mask_tmp = [mask_tmp(id_min:end,:); mask_tmp(1:id_min-1,:)];
                        if diff(mask_tmp(1:2,1)) > pi
                            % the azimuth is going counter clockwise, I want it clockwise with no
                            mask_tmp = [mask_tmp(1,:); mask_tmp(end:-1:2,:)];
                        end
                        flag_zenith = length(unique(sign(diff(mask_tmp(:,1))))) == 1; % the azimuth is monothonic crescent, meaning the mask cover zenith
                        az = mask_tmp(:,1) * 180 / pi;
                        mask = mask_tmp(:,2);
                        
                        % Case 2 patch crossing zenith
                        if flag_zenith
                            n_border = min(10, numel(az));
                            az = [az(end-n_border+1:end)-360; az; az(1:n_border)+360];
                            mask = [mask(end-n_border+1:end); mask; mask(1:n_border)];
                            % figure; plot(az_grid,interp1q(az, mask, az_grid));
                            % recompute mask_north
                            mask = interp1q(az, mask, az_grid) / pi * 180;
                            
                            for i = 1 : numel(az_grid)
                                out_mask(min(round(((90 - mask(i)) / el_step) + 1), size(out_mask,1)):-1:1, i) = true;
                            end
                        else
                            % Search the two minima
                            
                            for i = 1 : numel(az_grid)
                                id_ok = find(islocalmin(min(abs(az - az_grid(i)), 20)));
                                if numel(id_ok) == 2
                                    lim = sort(min(size(out_mask,1), round((pi/2 - mask(id_ok)) * 180/pi / el_step) + 1));
                                    out_mask(lim(1):lim(2), i) = true;
                                end
                            end
                        end
                    end
                    mask_tmp = mask_south;
                end
            end
            
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

        function [flag, weight] = isActive(this)
            % Get the status of activation of the constellation
            %
            % SYNTAX
            %   [flag weight] = isActive(this)
            flag = this.flag_enable;
            weight = this.weight;
        end

        function [weight] = getWeight(this)
            % Get the weight of the constellation
            %
            % SYNTAX
            %   [weight] = getWeight(this)
            weight = this.weight;
        end

        function setWeight(this, weight)
            % Set the weight of the constellation
            %
            % SYNTAX
            %   [weight] = setWeight(this)
            if not(isnumeric(weight)) || isempty(weight)
                Core.getLogger.addWarning(sprintf('Constellation weight for %s is not a valid number [0.0001 .. 100], using 1.0', this.SYS_NAME));
                weight = 1.0;
            elseif (weight < 0)
                Core.getLogger.addWarning(sprintf('Constellation weight for %s cannot be negative, using %.2f', this.SYS_NAME, weight));
                weight = abs(weight);
            end
            if (weight < 0.0001)
                Core.getLogger.addWarning(sprintf('Minimum weight for %s is 0.0001', this.SYS_NAME));
                weight = 0.0001;
            end
            if (weight > 100)
                Core.getLogger.addWarning(sprintf('Maximum weight for %s is 100', this.SYS_NAME));
                weight = 100;
            end
            this.weight = weight;
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
                this.setWeight(settings.getData(sprintf('%s_weight', this.SYS_NAME)));
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
                this.setWeight(settings.getWeight());
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
            str_cell = Ini_Manager.toIniString(sprintf('%s_weight', this.SYS_NAME), this.weight, str_cell);
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
    
    methods (Static)
        function [mask_north, mask_south] = generatePolarMask(i, r, lat, lon)
            % genarate polar mask all angle in radians
            %
            % SYNTAX:
            %    [mask1, mask2] = genratePolarMask(i,lat,lon)
            Z = r*sin(i);
            B = r*cos(i);
            alpha = (0:0.01:2*pi)';
            % mask north pole
            xyz_circle = [B*sin(alpha) B*cos(alpha) Z*ones(size(alpha))];
            [ox,oy,oz] = geod2cart(lat, lon, 0, GPS_SS.ELL_A, GPS_SS.ELL_F);
            xyz_circle(:,1) = xyz_circle(:,1) - ox;
            xyz_circle(:,2) = xyz_circle(:,2) - oy;
            xyz_circle(:,3) = xyz_circle(:,3) - oz;
            [circle_loc] = Coordinates.cart2local([ox,oy,oz], xyz_circle);
            hor_len = sqrt(circle_loc(:,1).^2 + circle_loc(:,2).^2);
            el = atan2(circle_loc(:,3),hor_len);
            el(circle_loc(:,3) < 0,:) = 0; % below horizon
            az = atan2(circle_loc(:,2),circle_loc(:,1)) - pi/2;
            mask_north = [az,el];
            if not(any(mask_north(:,2)))
                mask_north = [];
            end
            % mask south pole
            xyz_circle = [B*sin(alpha) B*cos(alpha) -Z*ones(size(alpha))];
            [ox,oy,oz] = geod2cart(lat, lon, 0, GPS_SS.ELL_A, GPS_SS.ELL_F);
            xyz_circle(:,1) = xyz_circle(:,1) - ox;
            xyz_circle(:,2) = xyz_circle(:,2) - oy;
            xyz_circle(:,3) = xyz_circle(:,3) - oz;
            [circle_loc] = Coordinates.cart2local([ox,oy,oz], xyz_circle);
            hor_len = sqrt(circle_loc(:,1).^2 + circle_loc(:,2).^2);
            el = atan2(circle_loc(:,3),hor_len);
            el(circle_loc(:,3) < 0,:) = 0; % below horizon
            az = atan2(circle_loc(:,2),circle_loc(:,1)) - pi/2;
            mask_south = [az,el]; 
            if not(any(mask_south(:,2)))
                mask_south = [];
            end
        end
    end
end
