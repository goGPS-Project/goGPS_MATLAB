%   CLASS Coordinates
% =========================================================================
%
% DESCRIPTION
%   Class to manage coordinates
%
% EXAMPLE
%   pos = Coordinates();
%
% CONSTRUCTOR SYNTAX
%   t = Coordinates.fromXYZ(XYZ);
%
% FOR A LIST OF CONSTANTs and METHODS use doc Coordinates
%
% COMMENTS
% The class stores arrays of time, not just a single element,
% it has been designed this way because MATLAB works best on arrays
%

%--- * --. --- --. .--. ... * ---------------------------------------------
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
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------


classdef Coordinates < Exportable & handle
    
    properties (Constant, GetAccess = public)
        ELL_A = GPS_SS.ELL_A;       % Ellipsoid semi-major axis
        ELL_E = GPS_SS.ELL_E;       % Ellipsoid eccentricity
        E2 = Coordinates.ELL_E^2;   % Square of the ellipsoidal eccentricity
        
        DEG2RAD = pi/180;           % Convert degree to radius
        RAD2DEG = 180/pi;           % Convert radius to degree
        
        VERSION = '2.0';            % New file version
        % 1.1 => adding s0_ip
        % 1.2 => adding observation rate
        % 1.3 => adding coo_type (fixed / non fixed)
        % 1.4 => adding master_name (name of the master station used as reference, one per epoch)
        % 1.5 => adding rate in header
        % 1.6 => new marker v3
        % 1.7 => adding reflectometry support
        % 1.8 => adding mean_snr
        % 2.0 => adding more reflectometry data (incompatible with version 1 - BREVA only)
        REF_TYPE = categorical({'ECEF','ENU','ENU_ROT'})
    end
    
    properties (SetAccess = public, GetAccess = public) % set permission have been changed from private to public (Giulio)
        name = ''                   % Name of the point
        name_v3 = ''                % 9 character long RINEX v3 name (it my or not be the same as name
        description = ''            % Point description
        
        time = GPS_Time             % Position time
        xyz = []                    % Coordinates are stored in meters in as cartesian XYZ ECEF [m]
        xyz_model = []              % Coordinates are stored in meters in as cartesian XYZ ECEF [m]
        v_xyz = []                  % Coordinates velocities XYZ ECEF  [m / year]
        discontinuity_id = []      % set discontinuities
        discontinuity               % discontinuity
        discontinuity_applied = false  %
        precision = 0.0001          % 3D limit [m] to check the equivalence among coordinates
        local_ref = struct('origin',[],'up_axis_rot',[])  % local coordinates system
        local_ref_type = categorical({'ECEF'})
        Cxx = []
        info = struct('n_epo', [], 'n_obs', [], 's0', [], 's0_ip', [], 'flag', [], 'fixing_ratio', [],'obs_used',[], 'rate', [], 'coo_type', [], 'mean_snr', [], 'master_name', uint64([])) % Additional info related to the coordinate in use
        reflectometry = struct('time', GPS_Time(), 'value', [], 'n_obs', [], 'std', []);
        rate = [];                  % coordinates rate: - default daily
        std_scaling_factor = 30;
    end
    
    % =========================================================================
    %% CONSTRUCTOR    
    % =========================================================================
    
    methods
        function this = Coordinates()
            % Constructor
            %
            % SYNTAX
            %   t = Coordinates();
        end
        
        function copyFrom(this, pos)
            % Copy from an object of the same type
            %
            % SYNTAX
            %   this.copyFrom(pos)
            try % legacy support
                this.name = pos.name;
            catch
                this.name = 'UNKN';
            end
            try % legacy support
                this.name_v3 = pos.name_v3;
            catch
                this.name_v3 = 'UNKN00WRD';
            end
            this.xyz = pos.xyz;
            this.Cxx = pos.Cxx;
            this.time = pos.time.getCopy;
            try % legacy support
                this.v_xyz = pos.v_xyz;
            catch
                this.v_xyz = [];
            end
            try % legacy support
                this.xyz_model = pos.xyz_model;
                this.discontinuity_id = pos.discontinuity_id;
                this.discontinuity = pos.discontinuity;
                this.discontinuity_applied = pos.discontinuity_applied;
            catch
                this.xyz_model = [];
                this.discontinuity_id = [];         % set discontinuities
                this.discontinuity;                  % discontinuity
                this.discontinuity_applied = false;  %
            end
            try % legacy support
                this.local_ref = pos.local_ref;
                this.local_ref_type = pos.local_ref_type;
            catch
                this.local_ref = struct('origin',[],'up_axis_rot',[]);  % local coordinates system
                this.local_ref_type = categorical({'ECEF'});
            end
            try % legacy support
                this.info = pos.info;
            catch
                this.info = struct('n_epo', [], 'n_obs', [], 's0', [], 's0_ip', [], 'flag', [], 'fixing_ratio', [],'obs_used',[], 'rate', [], 'coo_type', [], 'mean_snr', [], 'master_name', uint64([])); % Additional info related to the coordinate in use
            end
            try % legacy support
                this.rate = pos.rate;
            catch
                this.rate = [];
            end
            
            try % legacy support
                this.reflectometry = pos.reflectometry;
                this.reflectometry.time = pos.reflectometry.time.getCopy;
            catch
                this.reflectometry = struct('time', GPS_Time(), 'value', [], 'n_obs', [], 'std', [], 'p2n', []);
            end

        end
        
        function az = getAzimuth(coo1, coo2)
            % GETAZIMUTH calculates the azimuth between two coordinates.
            %   az = GETAZIMUTH(coo1, coo2) calculates the azimuth from coordinate
            %   coo1 to coordinate coo2.
            %
            %  INPUT:
            %       coo1 - Coordinate 1 object
            %       coo2 - Coordinate 2 object
            %
            %  OUTPUT:
            %       az - Azimuth from coo1 to coo2 in degrees (0 to 360).
            %
            %  SINTAX:
            %   az = coo2.getAzimuth(coo2);

            coo1 = coo1.getMedianPos;
            coo2 = coo2.getMedianPos;

            % Convert latitude and longitude to radians
            lat1 = (coo1.lat);
            lon1 = (coo1.lon);
            lat2 = (coo2.lat);
            lon2 = (coo2.lon);
        
            % Calculate the longitude difference
            dlon = lon2 - lon1;
        
            % Calculate azimuth using the formula: arctan2(sin(dlon)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon))
            az = atan2(sin(dlon)*cos(lat2), cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon));
        
            % Convert azimuth from radians to degrees
            az = rad2deg(az);
        
            % Adjust azimuth range from -180 to 180 to 0 to 360
            if az < 0
                az = az + 360;
            end
        end

        function [bsl, dh] = getBaseline(coo1, coo2)
            % OUTPUT:
            %    baseline      distance
            %    height_diff   height difference
            % 
            % DESCRIPTION: overload of '-' function
            %
            % SYNTAX
            %   [baseline, height_diff] = getBaseline(coo1, coo2)
            %
            
            loc_enu = coo1.getMedianPos.getLocal(coo2.getMedianPos);
            bsl = hypot(loc_enu(1), loc_enu(2));
            dh = loc_enu(3);
        end 
        
        function [bsl, dh] = getLongBaseline(coo_ref, coo_list)
            % Get long baseline (vincentyDistance)
            %
            % INPUT
            %   coo_ref        reference coordinate
            %   coo_list       list of relative coordinates
            %
            % OUTPUT
            %    baseline      distance
            %    height_diff   height difference
            %             
            % SYNTAX
            %   [baseline, height_diff] = getBaseline(coo_ref, coo_list)
            
            [lat, lon, h] = coo_ref.getMedianPos.getGeodetic();
            [lat_list, lon_list, h_list] = deal(nan(numel(coo_list),1));
            for c = 1:numel(coo_list)
                [lat_list(c), lon_list(c), h_list(c)] = coo_list(c).getMedianPos.getGeodetic();
            end
            bsl = coo_ref.vincentyDistance(lat/pi*180, lon/pi*180, lat_list/pi*180, lon_list/pi*180);
            dh = h-h_list;
        end 
        
        function copy = getCopy(this)
            % Get a copy of this
            %
            % SYNTAX
            %   copy = getCopy(this)
            copy(numel(this)) = Coordinates();
            for i = 1:numel(this)
                copy(i).copyFrom(this(i));
            end
        end
        
        function this = sort(this)
            % Sort the coordinates in ascending time
            %
            % SYNTEX
            %   this = sort(this)
            
            [id_sort] = this.time.sort();
            if not(issorted(id_sort))
                this.xyz = this.xyz(id_sort, :);
                
                if ~isempty(this.Cxx)
                    this.Cxx = this.Cxx(:, :, id_sort);
                end
                
                % XYZ model
                if not(isempty(this.xyz_model)) && numel(this.xyz_model(:)) == numel(this.xyz(:))
                    this.xyz_model = this.xyz_model(id_sort, :);
                else
                    this.xyz_model = [];
                end
                
                % Number of epocs
                if not(isempty(this.info.n_epo))
                    this.info.n_epo = this.info.n_epo(id_sort);
                end
                
                % Number of observations
                if not(isempty(this.info.n_obs))
                    this.info.n_obs = this.info.n_obs(id_sort);
                end
                
                % Sigma0 of the solution
                if not(isempty(this.info.s0))
                    this.info.s0 = this.info.s0(id_sort);
                end
                
                % Sigma0 of the initial (pre-processing) solution
                if not(isempty(this.info.s0_ip))
                    this.info.s0_ip(id_sort) = this.info.s0_ip(id_sort);
                end
                
                % Validity flag
                if not(isempty(this.info.flag))
                    this.info.flag(id_sort) = this.info.flag(id_sort);
                end
                
                % Fixing ratio
                if not(isempty(this.info.fixing_ratio))
                    this.info.fixing_ratio(id_sort) = this.info.fixing_ratio(id_sort);
                end
                
                % Encyclopedia
                if not(isempty(this.info.obs_used))
                    this.info.obs_used(id_sort) = this.info.obs_used(id_sort);
                end
                
                % Rate of the original observations
                if not(isempty(this.info.rate))
                    this.info.rate(id_sort) = this.info.rate(id_sort);
                end
                
                % Coordinate type (Bernese style: F: fixed / G: rover)
                if not(isempty(this.info.coo_type))
                    this.info.coo_type(id_sort) = char(this.info.coo_type(id_sort));
                end
                
                % Coordinate type (Bernese style: F: fixed / G: rover)
                if ~isfield(this.info, 'mean_snr')
                    % Old Coordinates does not have mean_snr
                    this.info.mean_snr = [];
                end
                if not(isempty(this.info.mean_snr))
                    this.info.mean_snr = this.info.mean_snr(id_sort, :);
                else
                    this.info.coo.mean_snr = [];
                end
                
                % Name of the master station saved as numeric namev3
                if not(isempty(this.info.master_name))
                    this.info.master_name(id_sort) = this.info.master_name(id_sort);
                end
                
                % Reflectometry have is own time to sort -----
                
                if isempty(this.reflectometry.time)
                    if numel(this.reflectometry.value) == this.time.length()
                        this.reflectometry.time = this.time.getCopy;
                    end
                end
                [id_sort] = this.reflectometry.time.sort();
                if not(isempty(this.reflectometry.value))
                    this.reflectometry.value = this.reflectometry.value(id_sort);
                end
                if not(isempty(this.reflectometry.n_obs))
                    this.reflectometry.n_obs = this.reflectometry.n_obs(id_sort);
                end
                if not(isempty(this.reflectometry.std))
                    this.reflectometry.std = this.reflectometry.std(id_sort);
                end
            end
        end
        
        function this = append(this, pos, flag_refl_only)
            % Append a Coordinates object into the this
            %
            % SYNTAX
            %   this = append(this, pos)

            if nargin < 3
                flag_refl_only = false;
            end
            flag_refl_only = flag_refl_only && (this.time.length > 0); % if there are no coordinates, add also coordinates            
            
            if not(isempty(pos)) && not(all(pos.isEmpty))
                if ~this.isEmpty && ~this.time.isEmpty && pos.time.isEmpty
                    error('Inconsistent positions');
                end
                
                % Get the times of the reflectometry observations
                refl_time_list = [this.getReflectometryTime; pos.getReflectometryTime];
                
                % Start exporting coordinates
                if not(flag_refl_only)
                    n_el_old = size(this.xyz, 1);
                    n_el_new = size(pos.xyz, 1);

                    time_list = [this.time; pos.time];
                    [sync_time, id_sync] = time_list.getSyncedTime();
                    if isempty(id_sync)
                        id_sync = nan(size(this.xyz, 1) + size(this.xyz, 2), 2);
                        id_sync(1:size(this.xyz, 1), 1) = 1 : size(this.xyz, 1);
                        id_sync(size(this.xyz, 1) + (1:size(pos.xyz, 1)), 2) = 1 : size(pos.xyz, 1);
                    else
                        this.time = sync_time;
                    end

                    this.xyz = mergeNumericData(id_sync, this.xyz, pos.xyz);
                    if ~this.isEmpty && size(this.Cxx,1) > 0
                        if isempty(pos.Cxx)
                            pos.Cxx = nan(3,3,pos.length, class(this.Cxx));
                        else
                            Cxx = nan(3,3,size(id_sync, 1), class(pos.Cxx));
                            Cxx(:, :, not(isnan(id_sync(:,1)))) = this.Cxx(:,:,noNaN(id_sync(:,1)));
                            Cxx(:, :, not(isnan(id_sync(:,2)))) = pos.Cxx(:,:,noNaN(id_sync(:,2)));

                            this.Cxx = Cxx;
                        end
                    elseif ~isempty(pos.Cxx)
                        this.Cxx = pos.Cxx(:,:,id_sync(:,2));
                    end

                    this.xyz_model = []; % empty the model on append, with new coordinates this must be updated

                    try
                        % Number of epocs
                        this.info.n_epo = mergeNumericData(id_sync, this.info.n_epo, pos.info.n_epo);

                        % Number of observations
                        this.info.n_obs = mergeNumericData(id_sync, this.info.n_obs, pos.info.n_obs);

                        % Sigma0 of the solution
                        this.info.s0 = mergeNumericData(id_sync, this.info.s0, pos.info.s0);

                        % Sigma0 of the initial (pre-processing) solution
                        this.info.s0_ip = mergeNumericData(id_sync, this.info.s0_ip, pos.info.s0_ip);

                        % Validity flag
                        this.info.flag = mergeNumericData(id_sync, nan2zero(this.info.flag), nan2zero(pos.info.flag));

                        % Fixing ratio
                        this.info.fixing_ratio = mergeNumericData(id_sync, this.info.fixing_ratio, pos.info.fixing_ratio);

                        % Encyclopedia
                        this.info.obs_used = mergeNumericData(id_sync, this.info.obs_used, pos.info.obs_used);

                        % Rate of the original observations
                        this.info.rate = mergeNumericData(id_sync, this.info.rate, pos.info.rate);

                        % Coordinate type (Bernese style: F: fixed / G: rover)
                        this.info.coo_type = char(mergeNumericData(id_sync, uint8(this.info.coo_type(:)), uint8(pos.info.coo_type(:))));

                        % Mean SNR values per constellation
                        if ~isfield(this.info, 'mean_snr')
                            % If the field does not exist create it empty (old Coordinates)
                            this.info.mean_snr = [];
                        end
                        m_size_old = size(this.info.mean_snr);
                        m_size_new = size(pos.info.mean_snr);
                        if m_size_old(1) < n_el_old
                            % If the old data is shorter than it should be
                            this.info.mean_snr = [this.info.mean_snr; zeros(n_el_old - m_size_old(1), m_size_old(2), 'single')];
                        end
                        n_col = max(size(this.info.mean_snr,2), size(pos.info.mean_snr,2));
                        this.info.mean_snr = [this.info.mean_snr zeros(n_el_old, n_col - m_size_old(2))];
                        pos.info.mean_snr = [pos.info.mean_snr zeros(n_el_new, n_col - m_size_new(2))];
                        this.info.mean_snr = single(mergeNumericData(id_sync, single(this.info.mean_snr), single(pos.info.mean_snr)));

                        % Name of the master station saved as numeric namev3
                        this.info.master_name = mergeNumericData(id_sync, this.info.master_name, pos.info.master_name);
                    catch ex
                        Core_Utils.printEx(ex);
                    end
                end
                try 
                    % Reflectometry
                    [sync_time, id_sync] = refl_time_list.getSyncedTime();

                    if isempty(id_sync)
                        id_sync = nan(size(this.reflectometry.value, 1) + size(this.reflectometry.value, 2), 2);
                        id_sync(1:size(this.reflectometry.value, 1), 1) = 1 : size(this.reflectometry.value, 1);
                        id_sync(size(this.reflectometry.value, 1) + (1:size(pos.reflectometry.value, 1)), 2) = 1 : size(pos.reflectometry.value, 1);
                    else
                        this.reflectometry.time = sync_time;
                    end

                    this.reflectometry.value = mergeNumericData(id_sync, this.reflectometry.value, pos.reflectometry.value);
                    this.reflectometry.n_obs = mergeNumericData(id_sync, this.reflectometry.n_obs, pos.reflectometry.n_obs);
                    this.reflectometry.std   = mergeNumericData(id_sync, this.reflectometry.std,   pos.reflectometry.std);
                catch ex
                    Core_Utils.printEx(ex);
                end
                
                this.sort();
                this.setRate(this.getRate()); % Update the rate if needed
            end
            
            function data_out = mergeNumericData(id_sync, data1, data2)
                % merge two numeric array according to id_sync
                if isempty(data1)
                    if size(data2,2) > 1
                        data1(1,noNaN(id_sync(:,1))) = nan;  % row array
                    else
                        data1(noNaN(id_sync(:,1)), 1) = nan; % column array
                    end
                end
                if isempty(data2)
                    if size(data1,2) > 1
                        data2(1, noNaN(id_sync(:,2))) = nan;  % row array
                    else
                        data2(noNaN(id_sync(:,2)),1) = nan; % column array
                    end
                end
                data_out = zeros(size(id_sync, 1), max(size(data1,2), size(data2,2)), class(data2));
                if any(noNaN(id_sync(:,1)))
                    data_out(not(isnan(id_sync(:,1))), :) = data1(noNaN(id_sync(:,1)),:);
                end
                if any(noNaN(id_sync(:,2)))
                    data_out(not(isnan(id_sync(:,2))), :) = data2(noNaN(id_sync(:,2)),:);
                end
            end
        end
        
        function remBadTimes(coo_list)
            for coo = coo_list(:)'
                pos_time = coo.time.getMatlabTime*86400;
                rate = median(diff(pos_time(:)));
                time_ko = mod(pos_time / rate - median(mod(pos_time / rate,1)), 1) > 1e-5;
                if any(time_ko)
                    Core.getLogger.addWarning(sprintf('Coordinate %s object contains %d/%d non regularly sampled values', coo.getNameV3, sum(time_ko), numel(time_ko)));
                    coo.rem(time_ko);
                end
            end
        end

        function keep(coo_list, time_start, time_stop)
            % Keep coordinates into the interval of timees
            %
            % SYNTAX
            %   this = keep(this, time_start, timee_stop)
            
            coo_list.rem(GPS_Time(0), time_start);
            coo_list.rem(time_stop, GPS_Time('2100-01-01'));
        end
        
        function rem(coo_list, idx, time_stop)
            % Remove coordinates into the idx
            %
            % SYNTAX
            %   this = rem(this, idx)
            %   this = rem(this, time_start, timee_stop)
            
            if nargin == 3
                time_start = idx;
            end
            for c = 1 : numel(coo_list)
                coo = coo_list(c);
                
                if nargin == 3
                    idx = coo.time > time_start & coo.time < time_stop;
                else
                    if any(idx)
                        if islogical(idx)
                            idx = find(idx);
                        end
                        time_start = coo.time.getEpoch(idx(1));
                        time_stop = coo.time.getEpoch(idx(end));
                    else
                       time_stop = [];
                    end
                end
                
                coo.xyz(idx,:) = [];
                
                if ~isempty(coo.Cxx)
                    coo.Cxx(:,:,idx) = [];
                end
                if ~isempty(coo.v_xyz)
                    coo.v_xyz(idx,:) = [];
                end
                if ~isempty(coo.xyz_model)
                    coo.xyz_model(idx,:) = [];
                end
                % INFO
                if ~isempty(coo.info.n_epo)
                    coo.info.n_epo(idx) = [];
                end
                if ~isempty(coo.info.n_obs)
                    coo.info.n_obs(idx) = [];
                end
                if ~isempty(coo.info.s0)
                    coo.info.s0(idx) = [];
                end
                if ~isempty(coo.info.s0_ip)
                    coo.info.s0_ip(idx) = [];
                end
                if ~isempty(coo.info.flag)
                    coo.info.flag(idx) = [];
                end
                if ~isempty(coo.info.fixing_ratio)
                    coo.info.fixing_ratio(idx) = [];
                end
                if ~isempty(coo.info.obs_used)
                    coo.info.obs_used(idx) = [];
                end
                if ~isempty(coo.info.rate)
                    coo.info.rate(idx) = [];
                end
                if ~isempty(coo.info.coo_type)
                    coo.info.coo_type(idx) = [];
                end
                if ~isfield(coo.info, 'mean_snr')
                    % Old Coordinates does not have mean_snr
                    coo.info.mean_snr = [];
                end
                if ~isempty(coo.info.mean_snr)
                    coo.info.mean_snr(idx,:) = [];
                end
                if ~isempty(coo.info.master_name)
                    if islogical(idx)
                        n_data = numel(idx);
                    else
                        n_data = max(idx);
                    end
                    if numel(coo.info.master_name) < n_data
                        tmp = coo.info.master_name;
                        mtmp = median(tmp); % median master
                        if isempty(mtmp)
                            mtmp = uint64(0);
                        end
                        coo.info.master_name = [tmp; mtmp * ones(n_data - numel(tmp),1,'uint64')];
                    end
                    
                    coo.info.master_name(idx) = [];
                end
                
                coo.time.remEpoch(idx);
                
                tr = coo.getReflectometryTime;
                if not(isempty(time_stop)) && not(isempty(tr)) && not(tr.isEmpty())
                    % they might not have the same time,
                    % recompute idx
                    idx = coo.reflectometry.time > time_start & coo.reflectometry.time < time_stop;
                
                    coo.reflectometry.time.remEpoch(idx);

                    if ~isfield(coo.reflectometry, 'value')
                        coo.reflectometry.value = [];
                        coo.reflectometry.n_obs = [];
                        coo.reflectometry.std = [];
                    end
                    if ~isempty(coo.reflectometry.value)
                        coo.reflectometry.value(idx) = [];
                    end
                    if ~isempty(coo.reflectometry.n_obs)
                        coo.reflectometry.n_obs(idx) = [];
                    end
                    if ~isempty(coo.reflectometry.std)
                        coo.reflectometry.std(idx) = [];
                    end
                else
                    % do nothing, reflectometry observations are independent
                end
            end
        end
        
        function remReflectometryEntry(coo_list, idx, time_stop)
            % Remove coordinates into the idx
            %
            % SYNTAX
            %   this = remReflectometryEntry(this, idx)
            %   this = remReflectometryEntry(this, time_start, timee_stop)
            if nargin == 3
                time_start = idx;
            end
            for c = 1 : numel(coo_list)
                coo = coo_list(c);
                coo.getReflectometryTime; 

                if nargin == 3
                    idx = coo.reflectometry.time > time_start & coo.reflectometry.time < time_stop;
                end

                if ~isempty(coo.reflectometry.value)
                    coo.reflectometry.value(idx) = [];
                end
                if ~isempty(coo.reflectometry.n_obs)
                    coo.reflectometry.n_obs(idx) = [];
                end
                if ~isempty(coo.reflectometry.std)
                    coo.reflectometry.std(idx) = [];
                end
                coo.reflectometry.time.remEpoch(idx);
            end
        end

        function check(this)
            % Raise an error and empty the coordinates if time and
            % coordinates have different dimensions
            if size(this.xyz) ~= this.time.length
                fprintf('WARNING Coordinates with unconsistent dimension found\n');
                this.time = GPS_Time;
                this.xyz = [];
                this.xyz_model = [];
                this.xyz = [];
                this.Cxx = [];
                this.discontinuity_id = [];
                this.discontinuity = [];
                this.discontinuity_applied = [];
                n_data = 0;
                this.info = struct('n_epo', zeros(n_data, 1, 'uint32'), 'n_obs', zeros(n_data, 1, 'uint32'), 's0', zeros(n_data, 1, 'single'), 's0_ip', zeros(n_data, 1, 'single'), 'flag', zeros(n_data, 1, 'uint8'), 'fixing_ratio', zeros(n_data, 1, 'single'), 'obs_used',  zeros(n_data, 1, 'single'), 'rate', zeros(n_data, 1, 'single'), 'coo_type', char(ones(n_data, 1, 'uint8')) * 'U', 'mean_snr', zeros(n_data, 1, 'single'), 'master_name', zeros(n_data, 1, 'uint64'));
                this.reflectometry = struct('time', GPS_Time(), 'value', zeros(n_data, 1, 'single'), 'n_obs', zeros(n_data, 1, 'uint16'), 'std', zeros(n_data, 1, 'single'));
            end
        end
    end
    
    % =========================================================================
    %%    GETTERS
    % =========================================================================
    
    methods
        function [coo, id_coo] = get(coo_list, name)
            % get a coordinate with a specific name
            %
            % SINTAX
            %    coo_list.get(name)
            
            coo = Coordinates();
            id_coo = [];
            if iscell(name)
                c = 0;
                for name = name(:)'
                    c = c+1;
                    [coo(c), id_coo(c)] = coo_list.get(name{1});                    
                end
            else
                if isnumeric(name)
                    % This is a numeric marker v3
                    name = Core_Utils.markerCode2MarkerName(name);
                end
                id_coo = [];
                i = 0;
                found = false;
                if numel(name) == 9 % this is marker name V3
                    while i < numel(coo_list)
                        i = i + 1;
                        if strcmp(name, coo_list(i).getNameV3())
                            found = true;
                            id_coo = [id_coo; i];
                        end
                    end
                else % this is a short name marker (assuming v2 (4 char)
                    name = pad(name(1:min(4,numel(name))), 4); % convert to 4 char
                    while i < numel(coo_list)
                        i = i + 1;
                        cur_name = coo_list(i).getName();
                        cur_name = pad(cur_name(1:min(4,numel(cur_name))), 4); % convert to 4 char
                        if strcmp(name, cur_name)
                            found = true;
                            id_coo = [id_coo; i];
                        end
                    end
                end
                if found
                    coo = coo_list(id_coo);
                end
            end
        end
        
        function [id_ref] = getIdRef(coo_list)
            % Get the id of the reference station
            % The reference station has coo_type = 'F';
            % Return empty if not found
            %
            % SYNTAX
            %   coo_list.getIdRef();
            %
            id_ref = [];
            for i = 1:numel(coo_list)
                if (median(coo_list(i).info.coo_type) == 'F') 
                    id_ref = [id_ref; i]; 
                end
            end
        end

        function [name, descr] = getName(coo_list)
            % Get name and description
            %
            % SYNTAX
            %   marker_name = this.getNameV3();
            
            name = '';
            descr = '';
            for p = 1 : numel(coo_list)
                name{p} = upper(coo_list(p).name);
                % In legacy coordinate the field description was not present
                try
                    descr{p} = coo_list(p).description;
                catch
                    % use name instead
                    descr{p} = name{p};
                end
                if isempty(name{p})
                    name{p} = 'UNKN'; % Unknown name
                end
                if isempty(descr{p})
                    descr{p} = name{p};
                end
            end
            
            if numel(name) == 1
                name = name{1};
                descr = descr{1};
            end
        end
        
        function [name_v3, id_num] = getNameV3(coo_list)
            % Get marker name in RINEX V3 format
            %
            % SYNTAX
            %   marker_name = this.getNameV3();
            for p = 1 : numel(coo_list)
                this = coo_list(p);
                if Core_Utils.isMarkerV3(this.name_v3) && ~strcmp(this.name_v3(7:9), 'WRD')
                    if nargout == 2
                        [this.name_v3, id_num{p}] = Core_Utils.getMarkerV3(this.name_v3, this);
                    end
                else
                    [this.name_v3, id_num{p}] = Core_Utils.getMarkerV3(this.getName, this);
                end
                name_v3{p} = upper(this.name_v3);
            end
            
            if numel(name_v3) == 1
                name_v3 = name_v3{1};
                if nargout == 2
                    id_num = id_num{1};
                end
            end
        end
        
        function time = getTime(coo_list)
            % Get the time of the coordinates
            %
            % SYNTAX
            %   time = coo_list.getTime()
            
            time = [coo_list.time]; 
            time = [time.getCopy];
        end
        
        function time = getReflectometryTime(this)
            % Get time of the reflectometry measurement
            %
            % SYNTAX
            %   time = getReflectometryTime(this)
            %
            try
                no_time = isempty(this.reflectometry.time) || (isnumeric(this.reflectometry.time) &&  isnumeric(this.reflectometry.time == 0));
            catch
                no_time = true;
            end
            time = GPS_Time();
            if no_time
                % this is a new time, empty reflectometry
                if ~isempty(this.reflectometry) && numel(this.reflectometry.value) <= this.time.length()
                    % If there are less reflector heights wrt times, cut the first times
                    if any(this.reflectometry.value) && numel(noZero(this.reflectometry.value)) < this.time.length()
                        Core.getLogger.addWarning(sprintf('Corrupted number of reflector heights for "%s"', this.getName));
                    end
                    this.reflectometry.time = this.time.getEpoch((1:numel(this.reflectometry.value)) + (this.time.length() - numel(this.reflectometry.value)));
                    time = this.reflectometry.time.getCopy();
                end
            else
                time = this.reflectometry.time.getCopy;
            end
        end

        function [id_sync, time_sync] = getSyncId(coo_list, n_hours, flag_now)
            % Get time synced among all the coordinates
            %
            % SYNTAX
            %    [id_sync, time_sync] = coo_list.getSyncId(n_hours, flag_now)
            n_coo = numel(coo_list);
            if nargin > 1 && not(isempty(n_hours)) && n_hours > 0
                flag_cut = true;
                time_all = coo_list.getTime().tail(n_hours*3600);
            else
                flag_cut = false;
                time_all = coo_list.getTime();
            end
            rate_half = median(coo_list.getRate, 'omitnan')/2;
            time_start = inf;
            if nargin == 3 && not(isempty(flag_now)) && flag_now
                time_stop = floor(now * 86400 / rate_half) * rate_half / 86400;
            else
                time_stop = -inf;
            end
            for t = 1 : n_coo
                if ~time_all(t).isEmpty()
                    time_start = min(time_start, round(time_all(t).first.getMatlabTime * 86400 / rate_half) * rate_half / 86400);
                    time_stop  = max(time_stop, round(time_all(t).last.getMatlabTime * 86400 / rate_half) * rate_half / 86400);
                end
            end
            full_time = GPS_Time((time_start : rate_half/43200 : time_stop)');
            time_all = [time_all(:); full_time];
            [time_sync, id_sync] = time_all.getSyncedTime();
            if flag_cut
                time_lim = GPS_Time(time_stop);
                time_lim = time_lim.addIntSeconds(-(n_hours-1)*3600);
                id_ok = time_sync >= time_lim;
                time_sync = time_sync.getEpoch(id_ok);
                id_sync = id_sync(id_ok, :);
            end
            id_sync = id_sync(:, 1:end-1);
        end
        
        function [err_status, time_sync] = getProcessingStatus(coo_list, n_hours, flag_now)
            % Get error status synced among all the coordinates
            %
            % SYNTAX
            %    [err_status, time_sync] = coo_list.getProcessingStatus(n_hours, flag_now)
            
            switch nargin
                case 1
                    [id_sync, time_sync] = coo_list.getSyncId();
                case 2
                    [id_sync, time_sync] = coo_list.getSyncId(n_hours);
                case 3
                    [id_sync, time_sync] = coo_list.getSyncId(n_hours, flag_now);
            end
            
            err_status = uint32(isnan(id_sync));
            for c = 1 : numel(coo_list)
                if numel(noNaN(id_sync(:,c))) > 1
                    id_diff = -max(noNaN(id_sync(:,c))) + noNaN(id_sync(:,c));
                    is_pp = coo_list(c).isPreProcessed;
                    err_status(~isnan(id_sync(:,c)),c) = ...
                        10 * serialize(uint32(isnan(coo_list(c).info.s0_ip(end + id_diff)))) .* ...
                        serialize(uint32(not(is_pp(end + id_diff)))) + ...
                        100 * serialize(uint32(isnan(coo_list(c).info.s0(end + id_diff)))) .* ...
                        serialize(uint32(coo_list(c).info.coo_type(end + id_diff) ~= 'F'));
                end
            end
        end

        function [qi, flag] = getQualityIndex(this)
            % Compute a quality index based on the median indicators
            %
            % SYNTAX
            %   [qi, flag] = coo.getQualityIndex();

            qi = [double(this.info.n_epo) .* double(this.info.rate) ...
                double(this.info.n_obs) .* double(this.info.rate) ...
                double(this.info.fixing_ratio) ...
                double(this.isNETOk | this.isPPPOk) ...
                double(this.info.s0_ip) ...
                double(this.info.s0*1e3)];
            % Scale n_obs by n_epo
            % buffers are screwing this parameter
            rate = max(1,this.time.getRate);
            qi(:,2) = double(qi(:,2)) .* double(rate) ./ double(qi(:,1));
            % Saturate n_epo with the number of epoch in a session
            qi(:,1) = min(rate, qi(:,1));
            % Normalize the quality index
            norm_thr = median(qi, 'omitnan');
            norm_thr(1) = norm_thr(1) * 0.9; % of epochs
            norm_thr(2) = norm_thr(2) * 0.8; % of epochs
            norm_thr(5) = norm_thr(5) * 0.7; % of s0 pre-processing
            norm_thr(6) = norm_thr(6) * 0.7; % of s0 network / PPP 
            qi = min(1, qi ./ norm_thr);
            qi = prod(qi,2);
            flag = qi > 0.4;
        end
        
        function [sigmas_status, time_sync] = getProcessingSigmas(coo_list, n_hours, flag_now)
            % Get sigmas synced among all the coordinates
            %
            % SYNTAX
            %    [err_status, time_sync] = coo_list.getProcessingSigmas(n_hours, flag_now)
            
            switch nargin
                case 1
                    [id_sync, time_sync] = coo_list.getSyncId();
                case 2
                    [id_sync, time_sync] = coo_list.getSyncId(n_hours);
                case 3
                    [id_sync, time_sync] = coo_list.getSyncId(n_hours, flag_now);
            end
            
            sigmas_status = single(isnan(id_sync));
            for c = 1 : numel(coo_list)
                if numel(noNaN(id_sync(:,c))) > 1
                    id_diff = -max(noNaN(id_sync(:,c))) + noNaN(id_sync(:,c));
                    is_ip = (coo_list(c).info.coo_type(end + id_diff) == 'F') & ...
                            ~(~isnan(coo_list(c).info.s0(end + id_diff))  &...
                            coo_list(c).info.s0(end + id_diff) < 0.5);
                    sigmas_status(~isnan(id_sync(:,c)),c) = ...
                        nan2zero(coo_list(c).info.s0_ip(end + id_diff) .* ...
                        single(is_ip)) + ...
                        nan2zero(coo_list(c).info.s0(end + id_diff) .* ...
                        single(~is_ip));
                end
            end

        end
        
        function [fixing_status, time_sync] = getProcessingFixingRatio(coo_list, n_hours, flag_now)
            % Get fixing ratio synced among all the coordinates
            %
            % SYNTAX
            %    [err_status, time_sync] = coo_list.getProcessingFixingRatio(n_hours, flag_now)
            
            switch nargin
                case 1
                    [id_sync, time_sync] = coo_list.getSyncId();
                case 2
                    [id_sync, time_sync] = coo_list.getSyncId(n_hours);
                case 3
                    [id_sync, time_sync] = coo_list.getSyncId(n_hours, flag_now);
            end
            
            fixing_status = single(isnan(id_sync));
            for c = 1 : numel(coo_list)
                if numel(noNaN(id_sync(:,c))) > 1
                    id_diff = -max(noNaN(id_sync(:,c))) + noNaN(id_sync(:,c));
                    fixing_status(~isnan(id_sync(:,c)),c) = ...
                        nan2zero(coo_list(c).info.fixing_ratio(end + id_diff));
                end
            end

        end
        
        function coo = getMedianPos(coo_list, type)
            % get the median of the coordinates
            %
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % SYNTAX
            %   coo = getMedianPos(this, type)
            
            if nargin == 1 || isempty(type)
                type = 'obs';
            end
            
            coo = Coordinates();
            for c = 1 : numel(coo_list)
                this = coo_list(c);
                try
                    if isempty(this.time) || this.time.isEmpty
                        tmp = Coordinates.fromXYZ(median(this.getXYZ(type), 1, 'omitnan'));
                    else
                        tmp = Coordinates.fromXYZ(median(this.getXYZ(type), 1, 'omitnan'), this.time.getCentralTime);
                    end
                    coo(c) = this.getCopy;
                    coo(c).time = tmp.time;
                    coo(c).xyz = tmp.xyz;
                    coo(c).xyz_model = coo(c).xyz;
                catch ex
                    Core.getLogger.addWarning('No data found in coordinate object');
                    coo(c) = Coordinates();
                end
            end
        end
        
        function coo = getMergedPos(coo)
            if numel(coo) > 1
                % getSyncIndexes
                times = [coo.getTime];
                [sync_time, id_sync] = times.getSyncedTime();

                % compute baricentric position
                [phi, lam, h_ellips] = deal(nan(numel(coo),1));
                for r = 1 : numel(coo)
                    [phi(r), lam(r), h_ellips(r)] = coo(r).getMedianPos.getGeodetic;
                end
                m_phi = mean(phi, 'omitnan');
                m_lam = mean(lam, 'omitnan');
                m_h_ellips = mean(h_ellips, 'omitnan');

                % Compute sync coordinates
                [phi, lam, h_ellips] = deal(nan(size(id_sync,1), size(id_sync,2)));
                for r = 1 : numel(coo)
                    [phi(noNaN(id_sync(:,r)), r), lam(noNaN(id_sync(:,r)), r), h_ellips(noNaN(id_sync(:,r)), r)] = coo(r).getGeodetic;
                    [p, l, h] = coo(r).getMedianPos.getGeodetic;
                    phi(:,r) = phi(:,r) - p;
                    lam(:,r) = lam(:,r) - l;
                    h_ellips(:,r) = h_ellips(:,r) - h;
                end
                
                phi = robAdj(nan2zero(phi), double(~isnan(phi))) + m_phi;
                lam = robAdj(nan2zero(lam), double(~isnan(lam))) + m_lam;
                h_ellips = robAdj(nan2zero(h_ellips), double(~isnan(h_ellips))) + m_h_ellips;
                
                coo = Coordinates.fromGeodetic(phi, lam, h_ellips, [], sync_time);
                coo.setName('MRGD');
            end
        end

        function coo = getElement(this, id_el)
            % get a copy of the coordinates relative to a specific id (or epoch)
            %
            % UNFINISHED FUNCITON:
            %   covariances and precisions are not extracted
            %
            % SYNTAX
            %   coo = getElement(this, id_el)
            
            try
                if isempty(this.time) || this.time.isEmpty
                    % Simplified coo
                    coo = Coordinates.fromXYZ(this.xyz(id_el, :));
                    if numel(this.xyz) == numel(this.xyz_model)
                        coo.xyz_model = this.xyz_model(id_el,:);
                    end
                    coo.setName(this.getName);
                else
                    % full coo
                    coo = this.getCopy;
                    full_el = 1:this.length();
                    full_el(id_el) = [];
                    coo.rem(full_el);
                end
            catch
                Core.getLogger.addWarning('Coordinates get Element out of bound');
                coo = Coordinates();
            end
        end

        function coo_list = subSample(coo_list, new_rate)
            % Subsample coordinates according to the nw rate [s]
            %
            % SYNTAX
            %   coo_list = coo_list.subSample(new_rate);

            % get all the times
            time_list = [coo_list.time];
            time_list = time_list.round(1e-4);
            % generate a common time at a new rate
            [time_sync, id_sync] = time_list.getSyncedTime();
            t_start = time_sync.first.getNominalTime(new_rate);
            t_stop = time_sync.last.getNominalTime(new_rate);
            ref_stop = t_stop.getRefTime(t_start.getMatlabTime);
            new_time = GPS_Time.fromRefTime(t_start.getMatlabTime, 0:new_rate:ref_stop);
            % sync all the times to the same rate
            time_list = [new_time, time_list];
            [time_sync, id_sync] = time_list.getSyncedTime();
            % filter coordinates
            for c = 1 : numel(coo_list)
                id_ok = noNaN(id_sync(~isnan(id_sync(:,1)),c+1));
                coo_list(c) = coo_list(c).getElement(id_ok);
            end
        end
        
        
        function [xyz, y, z, time] = getXYZ(this, type)
            % Get Coordinates as cartesian Earth Centered Earth Fixed Coordinates
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   xyz      = coordinates     [m]
            %
            % SYNTAX
            %   [xyz, time] = this.getXYZ()
            %   [x, y, z, time] = this.getXYZ()
            
            if nargin == 1 || isempty(type)
                type = 'obs';
            end
            
            if type(1) == 'm'
                xyz_tmp = this.xyz_model;
                if isempty(xyz_tmp)
                    xyz_tmp = zeros(size(this.xyz));
                end
            else
                xyz_tmp = this.xyz;
            end
            % Check if empty
            if isempty(xyz_tmp)
                xyz_tmp = [nan nan nan];
            end
            
            if nargout >= 3
                xyz = xyz_tmp(:, 1);
                y =   xyz_tmp(:, 2);
                z =   xyz_tmp(:, 3);
                if isempty(this.time)
                    this.time = GPS_Time();
                else
                    time = this.time.getCopy;
                end
            else
                xyz = xyz_tmp;
                if isempty(this.time)
                    this.time = GPS_Time();
                end
                y = this.time.getCopy;
            end
        end
        
        function [east, north, utm_zone, time] = getUTM(this, type)
            % Get Coordinates as UTM coordinates
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   east     = east Coordinates  [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone          [4char]
            %
            % SYNTAX
            %   [east, north, utm_zone, time] = this.getUTM();
            
            if nargin < 2 || isempty(type)
                type = 'obs';
            end
            
            [lat, lon] = this.getGeodetic(type);
            [east, north, utm_zone] = this.geod2utm(lat, lon);
            if nargout == 4
                time = this.time.getCopy;
            elseif nargout <=2 
                    east = [east; north];            
            end            
            if nargout == 2
                north = this.time.getCopy;
            end
        end
        
        function [east, north, up, utm_zone, time] = getENU(this, theta)
            % Get Coordinates as UTM ENUs coordinates
            %
            % INPUT
            %   theta   planar rotation angle [degree -360:360]
            %
            % OUTPUT
            %   east     = east Coordinates  [m]
            %   north    = north Coordinates [m]
            %   up       = up                [m]
            %   utm_zone = UTM zone          [4char]
            %
            % SYNTAX
            %   [east, north, up, utm_zone, time] = this.getENU(<theta>);
            %   [enu, utm_zone, time] = this.getENU(<theta>);
            
            [lat, lon, up] = this.getGeodetic();
            [east, north, utm_zone] = this.geod2utm(lat, lon);
            
            if nargin == 2 && ~isempty(theta)
                tmp = [east(:) north(:)] * [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
                east = tmp(:,1);
                north = tmp(:,2)';
            end
            if nargout <= 3
                east = [east(:) north(:) up(:)];
                north = utm_zone;
                if nargout == 3
                    up = this.time.getCopy;
                end
            end
            if nargout == 5
                time = this.time.getCopy;
            end
        end
        
        function h_ortho = getOrthometricHeight(this)
            % Get Coordinate orthometric height
            %
            % INPUT
            %
            % OUTPUT
            %   h_ortho  = orthometric height            [m]
            %
            % SYNTAX
            %   [h_ortho] = this.getOrthometricHeight()
            [lat, lon, h_ellips, h_ortho] = this.getGeodetic();
        end

        function h_ellips = getEllipsoidalHeight(this)
            % Get Coordinates ellipsoidal height
            %
            % INPUT
            %
            % OUTPUT
            %   h_ellips = ellipsoidal height            [m]
            %
            % SYNTAX
            %   [h_ellips] = this.getEllipsoidalHeight()
            [lat, lon, h_ellips, h_ortho] = this.getGeodetic();
        end

        function [lat, lon, h_ellips, h_ortho] = getGeodetic(this, type)
            % Get Coordinates as Geodetic coordinates
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   lat      = latitude                      [rad]
            %   lon      = longitude                     [rad]
            %   h_ellips = ellipsoidal height            [m]
            %   h_ortho  = orthometric height            [m]
            %
            % SYNTAX
            %   [lat, lon, h_ellips, h_ortho] = this.getGeodetic(type)
            
            if nargin < 2 || isempty(type)
                type = 'obs';
            end
            if nargout > 2
                [lat, lon, h_ellips] = Coordinates.cart2geod(this.getXYZ(type));
                if nargout == 4
                    h_ortho = h_ellips - this.getOrthometricCorrFromLatLon(lat, lon);
                end
            else
                [lat, lon] = Coordinates.cart2geod(this.xyz);
                if nargout <= 1
                    lat = [lat, lon];
                end
            end
        end
        
        function geod_lat = lat(this, type)
            % Get geodetic latitude
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   lat      = latitude                      [rad]
            %
            % SYNTAX
            %   [lat] = this.lat(<type>)
            
            if nargin < 2 || isempty(type)
                type = 'obs';
            end
            [geod_lat, geod_lon] = this.getGeodetic(type);
        end
        
        function geod_lon = lon(this, type)
            % Get geodetic latitude
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   lon      = longitude                     [rad]
            %
            % SYNTAX
            %   [lat] = this.lat(<type>)
            
            if nargin < 2 || isempty(type)
                type = 'obs';
            end
            [geod_lat, geod_lon] = this.getGeodetic(type);
        end
        
        
        function [lat_geoc, lon, h] = getGeocentric(this, type)
            % Get Coordinates as Geodetic coordinates
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   east     = east Coordinates [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone
            %
            % SYNTAX
            %   [lat, lon, h] = this.getGeocentric(type);
            if nargin < 2 || isempty(type)
                type = 'obs';
            end
            [lat_geoc, lon, h] = Coordinates.cart2geoc(this.getXYZ(type));
        end
        
        function ondu = getOrthometricCorrection(this, type)
            % Get Orthometric correction from the geoid loaded in Core
            %
            % INPUT
            %   type    it can be "obs" | "model"
            %           default "obs"
            %
            % OUTPUT
            %   ondu    = geoid ondulation [m]
            %
            % SYNTAX
            %   ondu = getOrthometricCorrection(this, type)
            
            if nargin < 2 || isempty(type)
                type = 'obs';
            end
            [lat, lon] = this.getGeodetic(type);
            ondu = this.getOrthometricCorrFromLatLon(lat, lon);
        end
        
        function [loc_enu, v_enu, id_ok] = getLocal(this, ref_pos, type)
            % Get Coordinates as Local coordinates with respect to ref_pos
            %
            % INPUT
            %   ref_pos     coo obj of the reference point, 
            %               if it is a series of coordinates takes the median
            %   type        it can be "obs" | "model"
            %               default "obs"
            %
            % OUTPUT
            %   loc_enu     xyz local coordinates
            %   v_enu       xyz local velocities
            %   id_ok       table of flags (outliers)
            %
            % SYNTAX
            %   [loc_enu, v_enu, id_ok] = this.getLocal(ref_pos, type)
            
            if nargin < 3 || isempty(type)
                type = 'obs';
            end
            this = this.getCopy;
            ref_pos = ref_pos.getCopy;
            tmp = [ref_pos; this(:)]; 
            tmp.setNewRef(1);
            ref_pos = tmp(1); 
            this = tmp(2:end);
            xyz_ref = ref_pos.getMedianPos.getXYZ;
            [xyz_this, time] = this.getXYZ(type);
            baseline = xyz_this - repmat(xyz_ref, size(xyz_this,1),1);
            [loc_enu, rot_mat] = Coordinates.cart2loca(xyz_ref, baseline);
            if nargout > 1
                if size(xyz_this, 1) > 3
                    if (nargout > 1)
                        v_enu = [0 0 0];
                        
                        for c = 1:3
                            tmp_time = time.getMatlabTime();
                            [~, id_ok(:,c), trend] = Coordinates.cooFilter([time.getMatlabTime loc_enu(:,c)], 0.8, 7);
                            id_ok(:,c) = not(id_ok(:,c)); % cooFilter returns KO epochs not OK
                            id_ok(:,c) = id_ok(:,c) & not(isnan(trend));
                            trend = Core_Utils.interp1LS(tmp_time(id_ok(:,c)), loc_enu(id_ok(:,c)), 1, tmp_time); % recompute trend with filtered data
                            v_enu(c) = (trend(end) - trend(1)) / (time.last.getMatlabTime - time.first.getMatlabTime) * 365; % m / year
                        end
                        if isempty(this.v_xyz)
                            this.v_xyz = v_enu * rot_mat;
                        end
                    end
                else
                    id_ok = true(size(xyz_ref,1));
                    v_enu = [nan nan nan];
                end
                
            end
        end
        
        function cov_xyz = getCovXYZ(this)
            % return variance covariance matrix in XYZ coordinates
            %
            % SYNTAX
            %   cov_xyz = this.getCovXYZ()
            
            cov_xyz = this.Cxx;
        end
        
        function cov_enu = getCovENU(this)
            % return variance covariance matrix in enu (local) coordinates
            %
            % SYNTAX
            %   cov_enu = this.getCovENU()
            
            [~, rot_mat] = Coordinates.cart2loca(this.getMedianPos.getXYZ, [0 0 0]);
            cov_enu = this.Cxx;
            for i = 1 :size(cov_enu,3)
                cov_enu(:,:,i) = rot_mat*cov_enu(:,:,i)*rot_mat';
            end
        end
        
        function std_xyz = getStdXYZ(this)
            % return std in XYZ coordinates
            %
            % SYNTAX
            %   std_xyz = this.getStdXYZ()
            
            std_xyz = nan(size(this.xyz,1),3);
            if ~isempty(this.Cxx)
                for i = 1 : size(this.Cxx,3)
                    std_xyz(i,:) = sqrt(diag(this.Cxx(:,:,i)));
                end
            end
        end
        
        function setStdXYZ(this, std_xyz)
            % return std in XYZ coordinates
            %
            % SYNTAX
            %   this.setStdXYZ(std_xyz)
            
            this.Cxx = zeros(3,3,size(this.xyz,1));
            if ~isempty(std_xyz)
                for i = 1 : size(this.Cxx,3)
                    this.Cxx(1,1,i) = real(std_xyz(i,1)^2);
                    this.Cxx(2,2,i) = real(std_xyz(i,2)^2);
                    this.Cxx(3,3,i) = real(std_xyz(i,3)^2);
                end
            end
        end
        
        function std_enu = getStdENU(this)
            % return vstd in ENU (local) coordinates
            %
            % SYNTAX
            %   cov_enu = this.getStdENU()
            
            [~, rot_mat] = Coordinates.cart2loca(this.getMedianPos.getXYZ, [0 0 0]);
            std_enu = nan(size(this.xyz,1),3);
            if ~isempty(this.Cxx)
                for i = 1 : size(this.Cxx,3)
                    std_enu(i,:) = real(sqrt(diag(rot_mat*this.Cxx(:,:,i)*rot_mat')));
                end
            end
        end
        
        function setStdENU(this, std_enu)
            % return std in ENU coordinates
            %
            % SYNTAX
            %   this.setStdENU(std_enu)
            
            [~, rot_mat] = Coordinates.cart2loca(this.getMedianPos.getXYZ, [0 0 0]);
            this.Cxx = zeros(3,3,size(this.xyz,1));
            if ~isempty(std_enu)
                for i = 1 : size(this.Cxx,3)
                    this.Cxx(1,1,i) = std_enu(i,1)^2;
                    this.Cxx(2,2,i) = std_enu(i,2)^2;
                    this.Cxx(3,3,i) = std_enu(i,3)^2;
                    this.Cxx(:,:,i) = rot_mat'*this.Cxx(:,:,i)*rot_mat; % ENU 2 XYZ
                end
            end
        end

        function status = isEmpty(coo_list)
            % Return the status of emptyness of the object
            %
            % SYNTAX
            %   status this.isEmpty();
            status = true(numel(coo_list), 1);
            for c = 1:numel(coo_list)
                status(c) = isempty(coo_list(c).xyz);
            end
        end
        
        function n_el = length(this)
            % Return the number of coordinates present in the object
            %
            % SYNTAX
            %   n_el = this.length();
            n_el = size(this.xyz, 1);
        end
        
        function is_pp = isPreProcessed(this, flag_list)
            % Find if the workspace have been pre-processed
            %
            % SYNTAX
            %    is_pp = this.isPreProcessed(<flag_list>)
            if nargin < 2
                flag_list = this.info.flag;
            end
            is_pp = bitand(flag_list, 1);
        end
        
        function is_ppp = isPPPOk(this, flag_list)
            % Find if the workspace have completed PPP
            %
            % SYNTAX
            %    is_ppp = this.isPPPOk()
            if nargin < 2
                flag_list = this.info.flag;
            end
            is_ppp = bitand(flag_list, 2);
        end
        
        function is_net = isNETOk(this, flag_list)
            % Find if the workspace have completed NET
            %
            % SYNTAX
            %    is_net = this.isNETOk()
            if nargin < 2
                flag_list = this.info.flag;
            end
            is_net = bitand(flag_list, 4);
        end
        
        function setPreProOk(this)
            % Set completed pre-processing
            %
            % SYNTAX
            %    is_pp = this.setPreProOk()
            this.info.flag = bitor(this.info.flag, 1);
        end
        
        function setPPPOk(this)
            % Set completed PPP
            %
            % SYNTAX
            %    is_pp = this.setNETOk()
            this.info.flag = bitor(this.info.flag, 2);
        end
        
        function setNETOk(this)
            % Set completed NET
            %
            % SYNTAX
            %    is_pp = this.setNETOk()
            this.info.flag = bitor(this.info.flag, 4);
        end
        function dist = sphDistanceTo(this, coo)
            % return the spherical distance on the ellipsoid between the object
            %
            % SYNTAX
            %     dist = this.sphDistanceTo(coo)
                        
            [lat1, lon1] = this.getGeodetic();
            [lat2, lon2] = coo.getGeodetic();
            r1 = sqrt(sum(this.getXYZ.^2));
            r2 = sqrt(sum(coo.getXYZ.^2));
            
            dist = Coordinates.sphericaldist(double(lat1) ./ pi * 180, double(lon1) ./ pi * 180, double(r1), double(lat2) ./ pi * 180, double(lon2) ./ pi * 180, double(r2));
        end
        
        function dist = ellDistanceTo(this, coo)
            % return the distance on the ellipsoid between the object
            % coordinates and the coordinate coo using vincenty's formula
            % (https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)
            %
            % SYNTAX
            %     dist = this.distanceTo(coo)
            
            % reduced latitude
            a = this.ELL_A;
            b = sqrt((1 - this.E2)*a^2) ;
            f = (a - b)/a;
            [lat1, lon1] = this.getGeodetic();
            [lat2, lon2] = coo.getGeodetic();
            lat1 = double(lat1);
            lat2 = double(lat2);
            U1 = atan((1-f)*tan(lat1));
            U2 = atan((1-f)*tan(lat2));
            L = lon2 - lon1;
            lambda_old = 10000;
            lambda = L;
            
            cos2_alpha = 0;
            while abs(lambda - lambda_old) > 1e-12
                sin_sigma = sqrt((cos(U2) * sin(lambda))^2 + (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda))^2);
                cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
                sigma = atan2(sin_sigma,cos_sigma);
                sin_alpha = cos(U1)*cos(U2)*sin(lambda)/sin_sigma;
                cos2_alpha = 1 - sin_alpha^2;
                cos_2sigmam = cos_sigma - 2*sin(U1)*sin(U2)/cos2_alpha;
                C = f/16*cos2_alpha*(4 + f*(4 -3*cos2_alpha));
                lambda_old = lambda;
                lambda = L + (1 - C)*f*sin_alpha*(sigma + C*sin_sigma*(cos_2sigmam + C*cos_sigma*(-1 + 2 * cos_2sigmam^2)));
            end
            u2 = cos2_alpha*(a^2 - b^2)/b^2;
            A = 1 + u2/16384*(4096 + u2*(-768 + u2*(320-175*u2)));
            B = u2/1024*(256 +u2*(-128 +u2*(74-47*u2)));
            Dsigma = B*sin_sigma*(cos_2sigmam + 1/4*B*(cos_sigma*(-1 + 2*cos_2sigmam^2) -B/6*cos_2sigmam*(-3 + 4 * sin_sigma^2)*(- 3 + 4 * cos_2sigmam^2)));
            dist = b*A*(sigma -Dsigma);
        end
        
        function [lid_ko, time, lid_ko_enu, trend_enu] = getBadSession(coo_list, coo_ref, n_obs)
            % Get outliers in East North Up coordinates
            %
            % SYNTAX
            %   [lid_ko lid_ko_enu, trend_enu] = getBadSession(coo_list, coo_ref, n_obs);
            %
            % SEE ALSO
            %   core.printKoSessions_experimental
            
            thr = 0.8;
            time = {};
            log = Core.getLogger();
            fh = figure('Visible', 'off'); Core_UI.beautifyFig(fh);
            for i = 1 : numel(coo_list)
                pos = coo_list(i);
                if ~pos.isEmpty
                    
                    if nargin == 1 || isempty(coo_ref)
                        enu_diff = pos.getLocal(pos.getMedianPos) * 1e3;
                        flag_time = true;
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            t = pos.time.getMatlabTime;
                            if numel(t) < size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more coordinates than times\n plotting only the positions with time'))
                                enu_diff = enu_diff(1:numel(t),:);
                            elseif numel(t) > size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more times than coordinates\n plotting only the first positions'))
                                t = t(1:size(enu_diff,1),:);
                            end
                        else
                            flag_time = false;
                            t = (1 : size(enu_diff, 1))';
                        end
                    elseif nargin > 1 % plot baseline
                        enu_diff = [];
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            [t_comm, idx_1, idx2] = intersect(round(coo_ref.time.getRefTime(pos.time.first.getMatlabTime)),round(pos.time.getRefTime(pos.time.first.getMatlabTime)));
                            t = pos.time.first.getMatlabTime + t_comm/86400;
                            enu_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz(idx2,:) - coo_ref.xyz(idx_1,:) )*1e3;
                            flag_time = true;
                        else
                            if numel(coo_ref.xyz) == numel(pos.xyz)
                                enu_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz - coo_ref.xyz)*1e3;
                                t = t(1:size(enu_diff,1),:);
                                flag_time = false;
                            else
                                log.addError(sprintf('No time in coordinates and number off coordinates in ref different from coordinate in the second receiver'))
                            end
                        end
                        enu_diff = bsxfun(@minus, enu_diff,median(enu_diff,1,'omitnan'));
                    end
                    
                    if nargin >= 2 && n_obs > 0
                        id_ok = (max(1, size(enu_diff,1) - n_obs + 1)) : size(enu_diff,1);
                        enu_diff = enu_diff(id_ok, :);
                        t = t(id_ok);
                    end
                    
                    e = enu_diff(1:numel(t),1);
                    [data, lid_ko_enu{i}(:,1), trend_enu{i}(:,1)] = Coordinates.cooFilter(e, 0.8, 7);
                    
                    n = enu_diff(:,2);
                    [data, lid_ko_enu{i}(:,2), trend_enu{i}(:,1)] = Coordinates.cooFilter(n, 0.8, 7);
                    
                    up = enu_diff(:,3);
                    [data, lid_ko_enu{i}(:,3), trend_enu{i}(:,1)] = Coordinates.cooFilter(up, 0.8, 7);
                    
                    lid_ko = isnan(e) | lid_ko_enu{i}(:,1) | isnan(n) | lid_ko_enu{i}(:,2) | isnan(up) | lid_ko_enu{i}(:,3);
                    time{i} = t;
                end
            end
            if numel(lid_ko) == 1
                lid_ko = lid_ko{1};
                lid_ko_enu = lid_ko_enu{1};
                trend_enu = trend_enu{1};
                time = time{1};
            end
        end
        
    end
    
    % =========================================================================
    %%    SETTERS
    % =========================================================================
    
    methods
        function setName(this, name, description)
            % Set the name of the coordinates
            %
            % SYNTAX
            %   this.setName(<name>,<description>)
            
            if not(isempty(name))
                name = upper(name);
                name_len = min(9,max(4, numel(name)));
                name = [name repmat(' ', 1, max(0, name_len - numel(name)))];
                this.name = name(1:name_len);
                if Core_Utils.isMarkerV3(name)
                    this.name_v3 = name;
                else
                    this.name_v3 = '';
                    this.name_v3 = this.getNameV3();
                end
            end
            if nargin == 3
                try
                    this.description = description;
                catch ex
                    % this try catch is here only for legacy support of old missing description field
                end
            end
            
        end
        
        function setTime(this, time)
            % Set the time of the coordinates
            %
            % SYNTAX
            %   this.setTime(time)
            
            if time.length > 3
                this.time = time.getNominalTime(time.getRate/2);
            else
                this.time = time.getCopy;
            end
            if this.time.length > size(this.xyz, 1)
                this.time.getEpoch(1 : size(this.xyz, 1));
                fprintf('The set coordinates time is larger than the number of positions stored\nCutting time\nDebug from Coordinates.setTime()');
            elseif this.time.length < size(this.xyz, 1)
                this.xyz = this.xyz(1 : this.time.length, :);
                this.time.getEpoch(1 : size(this.xyz, 1));
                fprintf('The set coordinates time is smaller than the number of positions stored\nCutting positions\nDebug from Coordinates.setTime()');
            end
            this.setRate(this.getRate);
        end
        
        function rate = getRate(coo_list)
            % Get the rate of the coordinates
            %
            % SYNTAX:
            %   rate = this.getRate
            
            rate = nan(numel(coo_list),1);
            for p = 1 : numel(coo_list)
                coo = coo_list(p);
                tmp = coo.rate;
                if isempty(tmp)
                    tmp = 1;
                end
                rate(p) = tmp;
                % Check if this is a good valid rate (consistent with 60% of the rates)
                flag_bad_rate = coo.time.length > 3 && sum(diff(sort(coo.time.getMatlabTime*86400)) == rate(p)) / (coo.time.length - 1) < 60;
                if flag_bad_rate || (isempty(coo.rate) ||  isnan(zero2nan(coo.rate))) && (not(isempty(coo.time)) && (coo.time.length > 2))
                    rate(p) = round(coo.time.getRate, 3);
                    coo.setRate(rate(p));
                end
            end
        end

        function sol_rate = getRateSolution(this)
            % Get the rate of the solution
            %
            % SYNTAX:
            %   sol_rate = this.getRateSolution
            
            sol_rate = this.getRate();
            if isempty(sol_rate) || isnan(zero2nan(sol_rate))
                sol_rate = this.rate;
                if isempty(sol_rate) || isnan(sol_rate)
                    sol_rate = median(diff(this.time.getMatlabTime*86400), 'omitnan');
                end
                if isempty(sol_rate) || isnan(sol_rate)
                    if numel(timestamp) < 2
                        sol_rate = Core.getState.sss_duration;
                    else
                        sol_rate = median(diff(timestamp * 86400), 'omitnan');
                    end
                end
            end
        end
        
        function setRate(this, rate)
            % Manually set the coordinate rate (do it carefully)
            %
            % SYNTAX
            %   this.setRate(rate);
            %
            % EXAMPLE
            %   this.setRate(Core.getState.sss_duration);
            this.rate = rate;
        end
        
        function setPosXYZ(this, xyz, y, z)
            % Set the Coordinates
            %
            % SYNTAX
            %   this.setPosXYZ(xyz)
            %   this.setPosXYZ(x, y, z)
            
            if nargin == 4
                xyz = [xyz(:) y(:) z(:)];
            end
            this.xyz = xyz;
        end
        
        function empty(this)
            % Empty the object
            %
            % SYNTAX
            %   this.empty();
            this.xyz = [];
            this.time = GPS_Time();
        end
        
        function addOffset(this, enu_offset, time_start, time_stop)
            % Add the offset (of the anntenna) to a set of coordinates
            %
            % SYNTAX
            %   this.addOffset(enu_offset, <time_start>, <time_stop>)
            if any(enu_offset)
                if nargin == 4
                    id_ok = (this.time >= time_start & this.time >= time_stop);
                else
                    id_ok = 1 : size(this.xyz,1);
                end
                
                xyz_offset = Coordinates.local2cart(this.getElement(id_ok).getMedianPos.xyz, enu_offset);
                this.xyz = this.xyz + repmat(xyz_offset, size(this.xyz,1), 1);
            end
        end
        
        function insertRHEntry(this, time, value, n_obs, std_val)
            % Insert in the reflectometry structure a value of reflector
            % height
            %
            % SYNTAX
            %   coo.insertRHEntry(time, value, n_obs, std_val)
            %   coo.insertRHEntry(reflectometry_struct)
            %
            
            if nargin == 2 && isstruct(time)
                value = time.value;
                n_obs = time.n_obs;
                std_val = time.std;
                time = time.time;
            end
            % Check coo reflectometry struct integrity
            t_ref = this.getReflectometryTime();
            
            if t_ref.length == 0                                
                this.reflectometry.time = time;
                this.reflectometry.value = single(value);
                this.reflectometry.n_obs = uint16(n_obs);
                this.reflectometry.std = single(std_val);
            else
                % find the closest time to the one of the inserted data

                for i = 1:time.length()
                    [t_dist, id_close] = min(abs(t_ref - time.getEpoch(i)));
                    if t_dist ~= 0 % this is not same time
                        % there is a new epoch!
                        this.reflectometry.time.append(time.getEpoch(i));
                        id_close = numel(this.reflectometry.value) + 1;
                    end
                    this.reflectometry.value(id_close,1) = single(value(i));
                    this.reflectometry.n_obs(id_close,1) = uint16(n_obs(i));
                    this.reflectometry.std(id_close,1) = single(std_val(i));
                end
            end
        end
    end
    
    % =========================================================================
    %%    STATIC CONSTRUCTOR
    % =========================================================================
    methods (Access = 'public', Static)
        function this = fromXYZ(xyz, y, z, time)
            % Set the Coordinates from XYZ coordinates
            %
            % SYNTAX
            %   this = Coordinates.fromXYZ(xyz)
            %   this = Coordinates.fromXYZ(x, y, z)
            
            this = Coordinates;
            if nargin > 2
                xyz = [xyz(:) y(:) z(:)];
                if nargin == 3
                    time = GPS_Time();
                end
            else
                if nargin == 2
                    time = y;
                else
                    time = GPS_Time();
                end
            end
            
            this.setPosXYZ(xyz);
            if ~time.isEmpty
                this.setTime(time);
            end
        end
        
        function this = fromCrdOut(name, crd_files, out_files)
            % Set the Coordinates from XYZ coordinates
            %
            % SYNTAX
            %   this = Coordinates.fromXYZ(xyz)
            %   this = Coordinates.fromXYZ(x, y, z)
            
            this = Coordinates;
            
            % read coordinates files
            logger = Logger.getInstance();
            
            fid = fopen(crd_files,'rt');
            
            if (fid < 0)
                logger.addError(['Failed to open ', crd_files]);
                time = [];
                xyz = [];
                enu = [];
                std_enu = [];
                std_3d = [];
            else
                txt = fread(fid,'*char')';
                logger.addMessage(['Reading ', crd_files]);
                fclose(fid);
                
                data = textscan(txt(392:end)',' %4d-%1d %4d-%03d %4d-%2d-%2d %1s %2d:%2d:%2d %2d:%2d:%2d %4s %6s %15f %15f %15f %15f %15f %15f %1s %8f %8f %8f %8f ');
                time = GPS_Time(datenum(double([data{:,[5:7]} data{:,[9:11]}])));
                xyz = [data{:,17:19}];
                enu = [data{:,20:22}];
                std_enu = [data{:,24:26}];
                std_3d = data{:,27};
                master = data{:,15};
            end
            fid = fopen(out_files,'rt');
            
            if (fid < 0)
                logger.addError(['Failed to open ', out_files]);
                n_epochs = [];
                n_l1 = [];
                enu = [];
                rms = [];
                rate = [];
            else
                txt = fread(fid,'*char')';
                logger.addMessage(['Reading ', out_files]);
                fclose(fid);
                
                data = textscan(txt(233:end)',' %4d-%1d %4d-%03d %4d-%2d-%2d %1s %2d:%2d:%2d %2d:%2d:%2d %8d %8d %8d %8d %8d %12f %8f');
                n_epochs = data{:,15};
                n_l1 = [];
                enu = [];
                rms = data{:,20};
                rate = data{:,21};
            end
            this.name = name;
            this.xyz = xyz;
            this.time = time;
            this.info.rate = rate;
            this.info.n_epo = n_epochs;
            this.info.master_name = master;
            this.info.s0_ip = rms;
            this.Cxx = zeros(3,3,size(xyz,1));
            for i = 1 : size(xyz,1)
                Cxxenu = zeros(3);
                Cxxenu(1,1) = std_enu(i,1);
                Cxxenu(2,2) = std_enu(i,2);
                Cxxenu(3,3) = std_enu(i,3);
                [lat, lon] = Coordinates.cart2geod(xyz(i,:));
                rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
                this.Cxx(:,:,i) = rot_mat*Cxxenu*rot_mat';
            end
            
        end
        
        function this = fromStringXYZ(xyz_string, time)
            % Set the Coordinates from XYZ coordinates (String)
            %
            % SYNTAX
            %   this = Coordinates.fromStringXYZ(xyz)
            
            this = Coordinates;
            xyz = sscanf(xyz_string, '%f%f%f')';
            if numel(xyz) ~= 3
                xyz = [0 0 0];
            end
            if nargin == 1
                time = GPS_Time();
            end
            this.setPosXYZ(xyz);
            this.setTime(time);
        end
        
        function this = fromGeodetic(lat, lon, h_ellips, h_ortho, time)
            % Set the Coordinates from Geodetic coordinates
            %
            % INPUT
            %   lat, lon [rad]
            %
            % SYNTAX
            %   this = Coordinates.fromGeodetic(phi, lam, h_ellips, [], time);
            %   this = Coordinates.fromGeodetic(phi, lam, [], h_ortho, time);
            
            this = Coordinates;
            if nargin >= 4 && ~isempty(h_ortho)
                h_ellips = h_ortho + this.getOrthometricCorrFromLatLon(lat, lon);
            end
            if nargin < 5
                time = GPS_Time();
            end            
            if ~((nargin > 2 && ~isempty(h_ellips)) || (nargin > 3 && ~isempty(h_ortho)))
                % No vertical coordinate => get it from SRTM
                h_ortho = Core_Utils.getElevation(lat, lon, true, {'etopo1'});
                h_ellips = h_ortho + this.getOrthometricCorrFromLatLon(lat, lon);                
            end
            N = GPS_SS.ELL_A ./ sqrt(1 - GPS_SS.ELL_E.^2 * sin(lat).^2);
            
            x = (N + h_ellips) .* cos(lon) .* cos(lat);
            y = (N + h_ellips) .* sin(lon) .* cos(lat);
            z = (N * (1 - GPS_SS.ELL_E.^2) + h_ellips) .* sin(lat);
            
            this = Coordinates.fromXYZ(x, y, z, time);
        end
        
        function this = fromUTM(easting, northing, utmzone, h_ortho, time)
            % Set the Coordinates from Geodetic coordinates
            %
            % INPUT
            %   lat, lon [rad]
            %
            % SYNTAX
            %   this = Coordinates.fromUTM(easting, northing, utmzone, h_ortho);

            [lat, lon] = Coordinates.utm2geod(easting, northing, utmzone);
            if nargin < 4 || isempty(h_ortho)
                h_ortho = zeros(size(easting));
            end
            if nargin == 5
                this = Coordinates.fromGeodetic(lat/180*pi, lon/180*pi, [], h_ortho, time);            
            else
                this = Coordinates.fromGeodetic(lat/180*pi, lon/180*pi, [], h_ortho);            
            end
        end
        
        function coo_list = fromFile(file_name)
            % Import from a coo file, mod file o mat file
            % File extension is used to discriminate the file type
            %
            % SYNTAX
            %   coo_list = Coordinates.fromFile(file_name);
            
            [folder, name, extension] = fileparts(file_name);
            coo_list = Coordinates();
            switch extension
                case {'.crd'}
                    coo_list = Coordinates.fromCrdFile(file_name);
                case {'.mod'}
                    coo_list = Coordinates.fromModFile(file_name);
                case {'.coo'}
                    coo_list = Coordinates.fromCooFile(file_name);
                case {'.mat'}
                    coo_list = Coordinates.fromMatFile(file_name);
                otherwise
                    Core.getLogger.addError(sprintf('Unknown file type "%s"', file_name));
            end
        end
        
        function coo_list = fromMatFile(file_name)
            % Import from a MATLAB file
            %
            % SYNTAX
            %   coo_list = Coordinates.fromMatFile(file_name);
            load(file_name, 'coo');
            coo_list = coo;
        end
        
        function coo_list = fromCooFile(file_name)
            % Import from a coo file XYZ and timestamp to a Coordinate object
            %
            % SYNTAX
            %   coo_list = Coordinates.fromCooFile(file_name);
            
            % cache for V3 names
            persistent cache4ch
            
            if not(iscell(file_name))
                file_name = {file_name};
            end
            
            coo_list = Coordinates;
            tmp_list = file_name;
            i = 0;
            file_name_list = {};
            for file_name = tmp_list
                file_name = file_name{1};
                file_info = dir(file_name);
                for f = 1 : numel(file_info)
                    i = i + 1;
                    file_name_list{i} = fullfile(file_info(f).folder, file_info(f).name);
                end
            end
            clear tmp_list
            f = 1;
            for file_name = file_name_list
                file_name = file_name{1};
                coo = Coordinates;
                
                if ~exist(file_name, 'file') == 2
                    % this is maybe a wild card
                end
                if exist(file_name, 'file') == 2
                    % Read and append
                    [txt, lim] = Core_Utils.readTextFile(file_name, 3);
                    if isempty(lim)
                        f = false;
                        timestamp = [];
                    else
                        % Verify the file version (it should be > 1.0):
                        id_ver = find(txt(lim(:,1) + 1) == 'F'); % +FileVersion
                        version_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), ['(?<=FileVersion[ ]*: )' coo.VERSION], 'once')));
                        file_ok = str2double(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), '(?<=FileVersion[ ]*: )[0-9].[0-9]', 'match', 'once')) >= 1;
                        
                        % Data should be present
                        timestamp = [];
                        data_start = size(lim, 1);
                        if file_ok
                            id_len_ok = find(lim(:,3) >= 9);
                            data_start = id_len_ok(find(txt(lim(id_len_ok,1) + 9) == 't') + 1); % +DataStart
                            id_len_ok = find(lim(:,3) >= 8);
                            data_stop = id_len_ok(find(txt(lim(id_len_ok,1) + 7) == 'd') -1); % +DataStop
                            if isempty(data_stop)
                                data_stop = size(lim, 1);
                            end
                            if isempty(data_start)
                                file_ok = false;
                            else
                                id_data = lim(data_start:data_stop,1);
                                % Read old timestamps
                                timestamp = datenum(txt(repmat(id_data, 1, 19) + repmat(0:18, numel(id_data), 1)), 'yyyy-mm-dd HH:MM:SS');
                            end
                        end
                        
                        if file_ok
                            % Point Name
                            id_line = find(txt(lim(1:data_start,1) + 1) == 'M'); % +MonitoringPoint
                            if isempty(id_line)
                                file_ok = false;
                            else
                                coo.name = regexp(txt(lim(id_line, 1):lim(id_line, 2)), '(?<=MonitoringPoint[ ]*: ).*', 'match', 'once');
                                coo.name_v3 = regexp(txt(lim(id_line, 1):lim(id_line, 2)), '(?<=NameV3[ ]*: ).*', 'match', 'once');
                                coo.description = regexp(txt(lim(id_line, 1):lim(id_line, 2)), '(?<=LongName[ ]*: ).*', 'match', 'once');
                            end
                            
                            % DataRate
                            id_line_rate = find(txt(lim(1:data_start-1,1) + 1) == 'D' & txt(lim(1:data_start-1,1) + 5) == 'R'); % +DateRate
                            if not(isempty(id_line_rate))
                                rate = str2double(regexp(txt(lim(id_line_rate,1) : lim(id_line_rate, 2)), '(?<=\:)[ 0-9\.]*', 'match', 'once'));
                                if isnan(rate)
                                    id_line_rate = []; % force recomputation
                                else
                                    coo.setRate(rate);
                                end
                            end
                            % DataType
                            id_line_start = find(txt(lim(1:data_start-1,1) + 1) == 'D' & txt(lim(1:data_start-1,1) + 5) == 'T'); % +DataType
                            id_line = id_line_start -1 + find(txt(lim(id_line_start:data_start-1,1) + 1) == '-');
                            col = str2num(txt(lim(id_line, 1) + repmat(2:3, numel(id_line),1))) + 1;
                            data_type = categorical();
                            for t = 1 : numel(col)
                                data_type(t) = categorical({txt((lim(id_line(t), 1) + 18) : lim(id_line(t), 2))});
                            end
                            
                            % Data column
                            % data_col = [2, 3, 4] + 1; % x, y, z
                            data_col = [col(data_type == categorical({'x'})), ...
                                col(data_type == categorical({'y'})), ...
                                col(data_type == categorical({'z'}))];
                        end
                        
                        
                        % Import data and time
                        if file_ok
                            coo.xyz = nan(data_stop - data_start + 1, 3);
                            n_data = data_stop - data_start + 1;
                            coo.info = struct('n_epo', zeros(n_data, 1, 'uint32'), 'n_obs', zeros(n_data, 1, 'uint32'), 's0', zeros(n_data, 1, 'single'), 's0_ip', zeros(n_data, 1, 'single'), 'flag', zeros(n_data, 1, 'uint8'), 'fixing_ratio', zeros(n_data, 1, 'single'), 'obs_used',  zeros(n_data, 1, 'single'), 'rate', zeros(n_data, 1, 'single'), 'coo_type', char(ones(n_data, 1, 'uint8')) * 'U', 'mean_snr', zeros(n_data, 1, 'single'), 'master_name', zeros(n_data, 1, 'uint64'));
                            coo.Cxx = zeros(3, 3, n_data);
                            
                            id_cov = [col(data_type == categorical({'Cxx'})), ...
                                col(data_type == categorical({'Cxy'})), ...
                                col(data_type == categorical({'Cxz'})), ...
                                col(data_type == categorical({'Cyy'})), ...
                                col(data_type == categorical({'Cyz'})), ...
                                col(data_type == categorical({'Czz'}))];
                            
                            id_n_epo = col(data_type == categorical({'nEpochs'}));
                            id_n_obs = col(data_type == categorical({'nObs'}));
                            id_s0_ip = col(data_type == categorical({'initialSigma0'}));
                            id_s0 = col(data_type == categorical({'sigma0'}));
                            id_fix = col(data_type == categorical({'fixingRatio'}));
                            id_rate = col(data_type == categorical({'obsRate'}));
                            id_ctype = col(data_type == categorical({'cooType'}));
                            id_master = col(data_type == categorical({'masterName'}));
                            for l = 0 : (data_stop - data_start)
                                [data_line] = regexp(txt(lim(data_start + l, 1) : lim(data_start + l, 2)), ';', 'split');
                                coo.xyz(l + 1, 1) = str2double(data_line{data_col(1)});
                                coo.xyz(l + 1, 2) = str2double(data_line{data_col(2)});
                                coo.xyz(l + 1, 3) = str2double(data_line{data_col(3)});
                                
                                if numel(id_cov) == 6
                                    tmp = str2double(data_line(id_cov));
                                    tmp = [tmp(1), tmp(2), tmp(3); ...
                                        tmp(2), tmp(4), tmp(5); ...
                                        tmp(3), tmp(5), tmp(6)]./1e6;
                                    if any(tmp(:))
                                        coo.Cxx(:,:,l + 1) = tmp;
                                    end
                                end
                                
                                if any(id_n_epo)
                                    coo.info.n_epo(l + 1) = uint32(str2double(data_line{id_n_epo}));
                                end
                                if any(id_n_obs)
                                    coo.info.n_obs(l + 1) = uint32(str2double(data_line{id_n_obs}));
                                end
                                if any(id_s0_ip)
                                    if id_s0_ip > numel(data_line)
                                        coo.info.s0_ip(l + 1) = nan;
                                    else
                                        coo.info.s0_ip(l + 1) = single(str2double(data_line{id_s0_ip}));
                                    end
                                end
                                if any(id_s0)
                                    if id_s0 > numel(data_line)
                                        coo.info.s0(l + 1) = nan;
                                    else
                                        coo.info.s0(l + 1) = single(str2double(data_line{id_s0}));
                                    end
                                end
                                if any(id_fix)
                                    if id_fix > numel(data_line)
                                        coo.info.fixing_ratio(l + 1) = nan;
                                    else
                                        coo.info.fixing_ratio(l + 1) = single(str2double(data_line{id_fix}));
                                    end
                                end
                                if any(id_rate)
                                    if id_rate > numel(data_line)
                                        coo.info.rate(l + 1) = nan;
                                    else
                                        coo.info.rate(l + 1) = single(str2double(data_line{id_rate}));
                                    end
                                end
                                if any(id_ctype)
                                    if id_ctype > numel(data_line)
                                        coo.info.coo_type(l + 1) = 'U';
                                    else
                                        coo.info.coo_type(l + 1) = char(data_line{id_ctype});
                                    end
                                end
                                if any(id_master)
                                    if id_master > numel(data_line)
                                        tmp = 'UNKN00WRD';
                                    else
                                        tmp = data_line{id_master};
                                    end
                                    if Core_Utils.isMarkerV3(tmp)
                                        % Just convert to numeric
                                        [~, coo.info.master_name(l + 1)] = Core_Utils.getMarkerV3(tmp);
                                    else
                                        % if it is not a marker V3 consider it a marker with 4 char
                                        tmp = pad(tmp(1:min(numel(tmp),4)), 4);
                                        id_4ch = Core_Utils.code4Char2Num(tmp);
                                        % search in cache
                                        id_cache = [];
                                        if not(isempty(cache4ch))
                                            id_cache = cache4ch(cache4ch(:,1) == uint64(id_4ch), 2);
                                        end
                                        if isempty(id_cache)
                                            [~, coo.info.master_name(l + 1)] = Core_Utils.getMarkerV3(tmp, coo);
                                            cache4ch = [cache4ch; uint64(id_4ch) uint64(coo.info.master_name(l + 1))];
                                        else
                                            coo.info.master_name(l + 1) = id_cache;
                                        end
                                    end
                                    
                                end
                            end
                            coo.time = GPS_Time(timestamp);
                            if isempty(id_line_rate)
                                coo.setRate(coo.getRate);
                            end
                        end
                        
                    end
                    
                    if isempty(coo.description)
                        coo.description = coo.name;
                    end
                    if not(Core_Utils.isMarkerV3(coo.name_v3))
                        coo.setName(coo.name); % this will set a V3 name;
                    end
                    
                    if not(version_ok)
                        % If the version is changed re-export the coordinates to update the file
                        log = Core.getLogger();
                        log.addMarkedMessage(sprintf('Update "%s" to the current Coordinates version %s', file_name, coo.VERSION));
                        coo.exportAsCoo(file_name);
                    end
                    coo_list(f) = coo;
                else
                    log = Core.getLogger();
                    log.addError(sprintf('%s cannot be imported', file_name));
                    coo_list(f) = Coordinates; % empty coordinates
                end
                f = f + 1;
            end
        end
        
        function coo_list = fromModFile(file_name)
            % Import from a model file XYZ, model and timestamp to a Coordinate object
            %
            % SYNTAX
            %   coo_list = Coordinates.fromCooFile(file_name);
            
            coo_list = Coordinates;
            if exist(file_name, 'file') == 2
                % Read and append
                [txt, lim] = Core_Utils.readTextFile(file_name, 3);
                if isempty(lim)
                    f = false;
                    timestamp = [];
                else
                    % Verify the file version (it should match 1.0):
                    version_ok = true;
                    file_ok = true;
                    
                    % Data should be present
                    timestamp = [];
                    data_start = size(lim, 1);
                    if file_ok
                        id_len_ok = find(lim(:,3) >= 9);
                        data_start = id_len_ok(find(txt(lim(id_len_ok,1) + 9) == 't') + 1); % +DataStart
                        id_len_ok = find(lim(:,3) >= 8);
                        data_stop = id_len_ok(find(txt(lim(id_len_ok,1) + 7) == 'd') -1); % +DataStop
                        if isempty(data_stop)
                            data_stop = size(lim, 1);
                        end
                        if isempty(data_start)
                            file_ok = false;
                        else
                            id_data = lim(data_start:data_stop,1);
                            % Read old timestamps
                            timestamp = datenum(txt(repmat(id_data, 1, 19) + repmat(0:18, numel(id_data), 1)), 'yyyy-mm-dd HH:MM:SS');
                        end
                    end
                    
                    if file_ok
                        % Point Name
                        id_line = find(txt(lim(1:data_start,1) + 1) == 'M'); % +MonitoringPoint
                        if isempty(id_line)
                            file_ok = false;
                        else
                            coo_list.name = regexp(txt(lim(id_line, 1):lim(id_line, 2)), '(?<=MonitoringPoint[ ]*: ).*', 'match', 'once');
                            coo_list.description = regexp(txt(lim(id_line, 1):lim(id_line, 2)), '(?<=LongName[ ]*: ).*', 'match', 'once');
                        end
                        
                        % DataRate
                        id_line_rate = find(txt(lim(1:data_start-1,1) + 1) == 'D' & txt(lim(1:data_start-1,1) + 5) == 'R'); % +DateRate
                        if not(isempty(id_line_rate))
                            rate = str2double(regexp(txt(lim(id_line_rate,1) : lim(id_line_rate, 2)), '(?<=\:)[ 0-9\.]*', 'match', 'once'));
                            if isnan(rate)
                                id_line_rate = []; % force recomputation
                            else
                                coo_list.setRate(rate);
                            end
                        end
                        % DataType
                        id_line_start = find(txt(lim(1:data_start-1,1) + 1) == 'D' & txt(lim(1:data_start-1,1) + 5) == 'T'); % +DataType
                        id_line = id_line_start -1 + find(txt(lim(id_line_start:data_start-1,1) + 1) == '-');
                        col = str2num(txt(lim(id_line, 1) + repmat(2:3, numel(id_line),1))) + 1;
                        data_type = categorical();
                        for t = 1 : numel(col)
                            data_type(t) = categorical({txt((lim(id_line(t), 1) + 18) : lim(id_line(t), 2))});
                        end
                        
                        % Data column
                        data_col = [col(data_type == categorical({'XMes'})), ...
                            col(data_type == categorical({'YMes'})), ...
                            col(data_type == categorical({'ZMes'}))];
                        data_model_col = [col(data_type == categorical({'XMod'})), ...
                            col(data_type == categorical({'YMod'})), ...
                            col(data_type == categorical({'ZMod'}))];
                    end
                    
                    % Import data and time
                    if file_ok
                        coo_list.xyz = nan(data_stop - data_start + 1, 3);
                        n_data = data_stop - data_start + 1;
                        coo_list.info = struct('n_epo', zeros(n_data, 1, 'uint32'), 'n_obs', zeros(n_data, 1, 'uint32'), 's0', zeros(n_data, 1, 'single'), 's0_ip', zeros(n_data, 1, 'single'), 'flag', zeros(n_data, 1, 'uint8'), 'fixing_ratio', zeros(n_data, 1, 'single'), 'obs_used',  zeros(n_data, 1, 'single'), 'rate', zeros(n_data, 1, 'single'), 'coo_type', char(ones(n_data, 1, 'uint8')) * 'U', 'mean_snr', zeros(n_data, 1, 'single'), 'master_name', zeros(n_data, 1, 'uint64'));
                        coo_list.Cxx = zeros(3, 3, n_data);
                        
                        for l = 0 : (data_stop - data_start)
                            [data_line] = regexp(txt(lim(data_start + l, 1) : lim(data_start + l, 2)), ';', 'split');

                            coo_list.xyz(l + 1, 1) = str2double(data_line{data_col(1)});
                            coo_list.xyz(l + 1, 2) = str2double(data_line{data_col(2)});
                            coo_list.xyz(l + 1, 3) = str2double(data_line{data_col(3)});
                            
                            coo_list.xyz_model(l + 1, 1) = str2double(data_line{data_model_col(1)});
                            coo_list.xyz_model(l + 1, 2) = str2double(data_line{data_model_col(2)});
                            coo_list.xyz_model(l + 1, 3) = str2double(data_line{data_model_col(3)});
                        end
                        coo_list.time = GPS_Time(timestamp);
                        if isempty(id_line_rate)
                            coo_list.setRate(coo_list.getRate);
                        end
                    end
                    
                    % Check description field
                    % use the old one for the file
                    coo.setName(coo.name); % this will set a V3 name;
                end 
            else
                log = Core.getLogger();
                log.addError(sprintf('%s cannot be imported', file_name));
            end
        end
        
        function coo_list = fromCrdFile(file_name)
            % Import from a breva crd file XYZ to a Coordinate object
            %
            % SYNTAX
            %   coo_list = Coordinates.fromCrdFile(file_name);
            
            coo_list = Coordinates();
            if ~isempty(file_name)
                [txt, lim] = Core_Utils.readTextFile(file_name);
                
                if isempty(txt)
                    Core.getLogger.addWarning(sprintf('File %s not found', file_name));
                    return
                end

                header_line = (txt(lim(:,1)) == '#');
                lim(header_line,:) = [];
                %initilaize array
                n_sta = size(lim,1);
                r= 0;
                for i = 1:n_sta
                    fprintf('Importing station %d/%d\n', i, n_sta);
                    line = txt(lim(i,1):lim(i,2));
                    parts = strsplit(line);
                    l = length(parts);
                    station_code = parts{1};
                    xyz = [str2double(parts{2}) str2double(parts{3}) str2double(parts{4})];
                    coo = Coordinates.fromXYZ(xyz(1), xyz(2), xyz(3), GPS_Time.now());
                    coo.setName(station_code);
                    flag_std = l >= 10 && ~isnan(str2double(parts{9}));
                    if l > 4
                        flag = str2double(parts{5});
                        if l > 7
                            coo.v_xyz = [str2double(parts{6}) str2double(parts{7}) str2double(parts{8})];
                            if flag_std
                                if l > 9
                                    std_pup = [str2double(parts{9}) str2double(parts{10})];
                                    if l > 11
                                        st_date = datenum([parts{11} ' ' parts{12}]);
                                        if not(isempty(st_date))
                                            coo.time = GPS_Time(st_date);
                                        end
                                        if l > 13
                                            en_date = datenum([parts{13} ' ' parts{14}]);
                                        end
                                    end
                                end
                            end
                        elseif l > 9
                            st_date = datenum([parts{9} ' ' parts{10}]);
                            if not(isempty(st_date))
                                coo.time = GPS_Time(st_date);
                            end
                            if l > 11
                                en_date = datenum([parts{11} ' ' parts{12}]);
                            end
                        end
                    end
                    if any(coo.xyz)
                        r = r + 1;
                        coo_list(r) = coo;
                    end
                end
            end
        end
        
        function coo_list = fromBernyCrdFile(file_name)
            % Import data from CRD file generated by Bernese processing
            % file_name can use wildcards e.g. filename = *.CRD
            %
            %  SYNTAX
            %   coo_list = Coordinates.fromBernyCrdFile(file_name);
            
            log = Core.getLogger.getInstance();
            
            file_list = dir(file_name);
            coo_list = Coordinates();
            
            for f = 1:numel(file_list)
                file_name = fullfile(file_list(f).folder, file_list(f).name);
                fid = fopen(file_name,'rt');
                
                if (fid < 0)
                    log.addError(['Failed to open ', file_name]);
                    coo_list(f) = Coordinates();
                else
                    txt = fread(fid,'*char')';
                    log.addMessage(['Reading ', file_name]);
                    fclose(fid);
                    
                    data = textscan(txt(392:end)',' %4d-%1d %4d-%03d %4d-%2d-%2d %1s %2d:%2d:%2d %2d:%2d:%2d %4s %6f %15f %15f %15f %15f %15f %15f %1s %8f %8f %8f %8f ');
                    time = GPS_Time(datenum(double([data{:,[5:7]} data{:,[9:11]}])));
                    rate = time.getRate;
                    if isempty(rate)
                        rate = 0;
                    end
                    time.addIntSeconds(rate/2);
                    xyz = [data{:,17:19}];
                    % enu = [data{:,20:22}];
                    std_enu = [data{:,24:26}];
                    % std_3d = data{:,27};
                    [dir_name, marker_name, file_ext] = fileparts(file_name);
                    
                    % Import into Coordinate object
                    coo_list(f) = Coordinates.fromXYZ(xyz, time);
                    coo_list(f).setName(marker_name);
                    coo_list(f).setStdENU(nan2zero(std_enu * 1e-3));
                    
                    coo_list(f).info.s0_ip = [data{:,16}];
                    
                    % Import master name code
                    for t = 1:time.length
                        name = data{1,15}{t};
                        if numel(name) < 9
                            name_tmp = '____00WRD'; % simulate NameV3
                            name_tmp(1:min(9,numel(name))) = name(1:min(9,numel(name)));
                            name = name_tmp;
                        end
                        coo_list(f).info.master_name(t) = Core_Utils.markerName2MarkerCode(name);
                    end
                end
            end
        end
    end
    
        
    % =========================================================================
    %%    INFOs
    % =========================================================================
    
    methods (Access = 'public')
        function str_description = toString(coo_list, flag_reverse)
            if nargin < 2 || isempty(flag_reverse)
                flag_reverse = true;
            end
            i = 0;
            for coo = coo_list
                i = i + 1;
                [lat, lon, h_ellips, h_ortho] = coo.getMedianPos.getGeodetic(); 
                lat = lat.*180/pi; lon = lon .*180/pi;
                [E,N,Z] = coo.getUTM();
                
                str = {sprintf('Master station: %s', coo.name), ...
                    iif(length(coo.description) <= 4, sprintf('%s', coo.getNameV3), sprintf('%s', coo.description)), ...
                    sprintf('Lat: %.6f deg', lat), ...
                    sprintf('Lon: %.6f deg', lon), ...
                    sprintf('El:  %.3f m', h_ellips), ...
                    sprintf('El:  %.3f m (orthometric)', h_ortho), ...
                    sprintf('UTM:  %.3f m E, %.3f m N, %s', E, N, Z)};

                if flag_reverse
                    loc_info = Core_Utils.getLocationInfo(lat, lon);
                    if ~isempty(loc_info)
                        str = [str, {'--------------------------------', ...
                            sprintf('%s', loc_info.display_name)}];
                    end
                end
                str_description{i} = '';
                for str_line = str
                    str_description{i} = sprintf('%s%s\n',str_description{i}, str_line{1});
                end
            end
            if nargout == 0
                for i = 1:numel(coo_list)
                    fprintf('%s\n', str_description{i});
                end               
            end
        end

        function printSessagesimalCoordinates(coo_list)
            if not(isempty(coo_list))
                fprintf('\n------------------------------------------------------------------\n');
                fprintf('MARKER NAME  |   Lat [deg] |   Lon [deg] | h_ellips [m] | h_ortho [m] \n')
                fprintf('------------------------------------------------------------------\n');
            else
                fprintf('No data found');
            end
            for coo = coo_list
                [lat, lon, h_ellips, h_ortho] = coo.getGeodetic();
                n_ch = 13;
                name = coo(1).name;
                if numel(name) < n_ch
                    name = [name char(32*ones(1, n_ch-numel(name), 'uint8'))];
                else
                    name = name(1:13);
                end
                lat = lat/pi*180;
                lon = lon/pi*180;
                lat_sess = nan(numel(lat),3);
                lat_sess(:,1) = fix(lat);
                lat_sess(:,2) = fix((lat-lat_sess(:,1))*60);
                lat_sess(:,3) = round((lat-lat_sess(:,1)-lat_sess(:,2)/60)*60*60);
                lon_sess = nan(numel(lat),3);
                lon_sess(:,1) = fix(lon);
                lon_sess(:,2) = fix((lon-lon_sess(:,1))*60);
                lon_sess(:,3) = round((lon-lon_sess(:,1)-lon_sess(:,2)/60)*60*60);
                for i = 1:numel(lat)
                    fprintf('%s| %3d %2d'' %2d" | %3d %2d'' %2d" |%13.2f |%12.2f\n', name, lat_sess(i,1), lat_sess(i,2), lat_sess(i,3), lon_sess(i,1), lon_sess(i,2), lon_sess(i,3), h_ellips, h_ortho);
                end
            end
        end

        function printTable(coo_list, flag_info)
            if nargin == 1 || isempty(flag_info)
                flag_info = false;
            end
            if not(isempty(coo_list))
                if flag_info
                    fprintf('\n---------------------------------------------------------------------------------------------------------------------------------------------\n');
                    fprintf('MARKER NAME  |  Lat [deg] |  Lon [deg] | h_ellips [m] | h_ortho [m] | Address\n')
                    fprintf('---------------------------------------------------------------------------------------------------------------------------------------------\n');
                else
                    fprintf('\n------------------------------------------------------------------\n');
                    fprintf('MARKER NAME  |  Lat [deg] |  Lon [deg] | h_ellips [m] | h_ortho [m] \n')
                    fprintf('------------------------------------------------------------------\n');
                end
            else
                fprintf('No data found');
            end
            
            for coo = coo_list
                [lat, lon, h_ellips, h_ortho] = coo.getMedianPos.getGeodetic();
                n_ch = 13;
                name = coo(1).name;
                if numel(name) < n_ch
                    name = [name char(32*ones(1, n_ch-numel(name), 'uint8'))];
                else
                    name = name(1:13);
                end
                if flag_info && not(isnan(lat))
                    loc_info = Core_Utils.getLocationInfo(lat/pi*180, lon/pi*180);
                    fprintf('%s|%11.6f |%11.6f |%13.2f |%12.2f | %s\n', name, lat/pi*180, lon/pi*180, h_ellips, h_ortho, loc_info.display_name);
                else
                    fprintf('%s|%11.6f |%11.6f |%13.2f |%12.2f\n', name, lat/pi*180, lon/pi*180, h_ellips, h_ortho);
                    
                end
            end
            fprintf('\n');
        end
        
        function dark_sky_link = openDarkSky(coo, date)
            % Get darkSky link (and open it in a browser)
            %
            % INPUT
            %   coo     single coordinate object (not an array)
            %   date    (optional) if empty show the page of the forecast
            %                      if not provided as parameter show the
            %                      last day contained into the coo object
            %
            % SYNTAX
            %   dark_sky_link = openDarkSky(coo, <date>)
            [lat, lon] = coo.getMedianPos.getGeodetic();
            if nargin == 1
                if not(isempty(coo.time))
                    date = coo.time.last;
                else
                    date = [];
                end
            end
            if not(isempty(date))
                dark_sky_link = sprintf('https://darksky.net/details/%.4f,%.4f/%s/ca12/en', lat .* (180 / pi), lon .* (180 / pi), date.toString('yyyy-mm-dd'));
            else
                dark_sky_link = sprintf('https://darksky.net/forecast/%.4f,%.4f/ca12/en', lat .* (180 / pi), lon .* (180 / pi));
            end
            web(dark_sky_link);
        end
        
        function refl_link = openReflectionZoneMapping(coo_list, refl_height, az_lim);
            % Open the webpage of gnss-reflections to show
            % reflection the reflection zone mapping
            %
            % INPUT
            %   coo           single coordinate object (not an array)
            %   refl_height   (optional) antenna height above reflecting surface [m]
            %                 default = 1.5 m
            %
            % SYNTAX
            %   refl_link = openReflectionZoneMapping(coo, <ant_height>)
            
            if nargin < 3
                az_lim = [0 360];
            end
            
            for r = 1 : numel(coo_list)
                [lat, lon,h_ellips] = coo_list(r).getMedianPos.getGeodetic();
                if nargin == 1
                    refl_height = 1.5;
                end
                refl_link{r} = sprintf('http://gnss-reflections.org/rzones?station=%s&lat=%.4f&lon=%.4f&height=%.3f&msl=off&RH=%.3g&eang=3&azim1=%d&azim2=%d', ...
                    coo_list(r).getName(), ...
                    lat .* (180 / pi), lon .* (180 / pi), ...
                    h_ellips, ...
                    refl_height, az_lim(1), az_lim(2));
                web(refl_link{r});
            end
        end   
        
        function geoid_link = openGeoidPage(coo)
            % Open the webpage of gnss-reflections to show 
            % a google map and geodetic height of the point
            %
            % INPUT
            %   coo           single coordinate object (not an array)
            %
            % SYNTAX
            %   refl_link = openReflectionZoneMapping(coo, <ant_height>)
            [lat, lon,h_ellips] = coo.getMedianPos.getGeodetic();
            if nargin == 1
                ant_height = 1.5;
            end
            geoid_link = sprintf('http://gnss-reflections.org/geoid?station=%s&lat=%.4f&lon=%.4f&height=%.3f', ...
                coo.getName(), ...
                lat .* (180 / pi), lon .* (180 / pi), ...
                h_ellips);
            web(geoid_link);
        end
    end
    
    % =========================================================================
    %%    SHOW
    % =========================================================================
    
    methods (Access = 'public')
        function fh = showBaselineENU(coo_list, id_ref, id_list_sta, n_pts, rot_angle_cw)
            % Plot East North Up Baseline with respect to the coordinate with id = id_ref
            %
            % SYNTAX
            %   coo_list.showBaselineENU(id_ref, id_ref, id_list_sta, n_pts, rot_angle_cw);
            if (nargin < 2) || isempty(id_ref)
                id_ref = 1:numel(coo_list);
            else
                if ischar(id_ref)
                    [~, id_ref] = coo_list.get(id_ref);
                end
            end
            if (nargin < 3) || isempty(id_list_sta)
                id_list_sta = 1:numel(coo_list);
            else
                if ischar(id_list_sta)
                    [~, id_list_sta] = coo_list.get(id_list_sta);
                elseif iscell(id_list_sta)
                    tmp = id_list_sta;
                    id_list_sta = [];
                    for i = 1:numel(tmp)
                        [~, tmp2] = coo_list.get(tmp{i});
                        if ~isempty(tmp2)
                            id_list_sta = [id_list_sta(:); tmp2];
                        end
                    end
                end
            end
            if nargin < 5 || isempty(rot_angle_cw)
                rot_angle_cw = 0;
            end
            if (nargin >= 4) && ~isempty(n_pts)
                if numel(id_ref) == 1
                    fh = coo_list(setdiff(id_list_sta, id_ref)).showCoordinatesENU(coo_list(id_ref), n_pts, rot_angle_cw);
                else
                    fh = coo_list(id_list_sta).showCoordinatesENU(coo_list(id_ref), n_pts, rot_angle_cw);
                end
            else
                if numel(id_ref) == 1
                    fh = coo_list(setdiff(id_list_sta, id_ref)).showCoordinatesENU(coo_list(id_ref), 0, rot_angle_cw);
                else
                    fh = coo_list(id_list_sta).showCoordinatesENU(coo_list(id_ref), 0, rot_angle_cw);
                end
            end
        end

        function fh = showBaselinePlanarUp(coo_list, id_ref)
            % Plot Planar and Up Baseline with respect to the coordinate with id = id_ref
            %
            % SYNTAX
            %   coo_list.showBaselinePlanarUp(id_ref);
            if (nargin < 2) || isempty(id_ref)
                id_ref = 1:numel(coo_list);
            else
                if ischar(id_ref)
                    [~, id_ref] = coo_list.get(id_ref);
                end
            end
            %if numel(id_ref) > 1
            fh = coo_list.showCoordinatesPlanarUp(coo_list(id_ref));
            %else
            %    fh = coo_list.showCoordinatesPlanarUp(coo_list(id_ref));
            %end
        end

        function fh = showPositionENU(coo_list)
            % Plot East North Up coordinates
            %
            % SYNTAX
            %   this.showPositionENU(coo_list);
            fh = showCoordinatesENU(coo_list);
        end
        
        function fh_list = showNData(coo_list)
            % Simple display of the number of data /obs used to compute the solution
            %
            % SYNTAX
            %   fh_list = coo_list.showNData();
            fh_list = [];
            for coo = coo_list(:)'
                if ~isempty(coo.info.n_obs)
                    if not(isempty(coo.name))
                        fig_name = sprintf('%s #data', coo.name);
                    else
                        fig_name = sprintf('#data');
                    end
                    fh = figure('Visible', 'off');  Core_UI.beautifyFig(fh);
                    fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                    fh.UserData = struct('fig_name', fig_name);

                    setAxis(fh);
                    plotSep(coo.time.getMatlabTime, coo.info.n_obs, '.-', 'MarkerSize', 10, 'LineWidth', 2);
                    ylabel('n obs');
                    setAxis(fh); yyaxis right
                    plotSep(coo.time.getMatlabTime, coo.info.n_epo, '.-', 'MarkerSize', 10, 'LineWidth', 2);
                    ylabel('n epochs');
                    title(sprintf('Number of data for %s\\fontsize{5} \n', coo.getName()));
                    xlim([coo.time.first.getMatlabTime coo.time.last.getMatlabTime]);
                    setTimeTicks();
                    
                    fh_list = [fh_list fh];
                    Core_UI.beautifyFig(fh);
                    Core_UI.addBeautifyMenu(fh);
                    fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                end
            end
        end
        
        function fh_list = showMeanSNR(coo_list)
            % Simple display of the meanSNR per receiver
            %
            % SYNTAX
            %   fh_list = coo_list.showMeanSNR();
            fh_list = [];
            for coo = coo_list(:)'
                if ~isempty(coo.info.n_obs)
                    if not(isempty(coo.name))
                        fig_name = sprintf('%s mean_snr', coo.name);
                    else
                        fig_name = sprintf('mean_snr');
                    end
                    fh = figure('Visible', 'off');  Core_UI.beautifyFig(fh);
                    fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                    fh.UserData = struct('fig_name', fig_name);

                    setAxis(fh);
                    if ~isfield(coo.info, 'mean_snr')
                        % Old Coordinates does not have mean_snr
                        coo.info.mean_snr = [];
                    end
                    plotSep(coo.time.getMatlabTime, zero2nan(coo.info.mean_snr), '.-', 'MarkerSize', 10, 'LineWidth', 2);
                    name = Constellation_Collector.SYS_NAME;
                    legend(name(1:size(coo.info.mean_snr,2)));
                    ylim([min(min(coo.info.mean_snr(:)), 15), max(max(coo.info.mean_snr(:)), 60)]);
                    ylabel('Mean SNR');
                    title(sprintf('Mean session SNR for %s\\fontsize{5} \n', coo.getName()));
                    xlim([coo.time.first.getMatlabTime coo.time.last.getMatlabTime]);
                    setTimeTicks();
                    
                    fh_list = [fh_list fh];
                    Core_UI.beautifyFig(fh);
                    Core_UI.addBeautifyMenu(fh);
                    fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                end
            end
        end
        
        function fh_list = showSigmas(coo_list, flag_prepro)
            % Simple display of the sigmas of the system used to compute the solution
            %
            % SYNTAX
            %   fh_list = coo_list.showSigmas();
            if nargin == 1
                flag_prepro = false;
            end
            fh_list = [];
            coo_list = coo_list(not(coo_list.isEmpty));
            
            yl = [0 -inf];
            xl = [inf -inf];
            
            str_title = 'Sigmas';
            if not(isempty(coo_list))
                fh = figure('Visible', 'off');  Core_UI.beautifyFig(fh);
            end
            name_list = {};
            for coo = coo_list(:)'
                if not(isempty(coo.name))
                    name_list = [name_list {strrep(coo.name, '_', ' ')}];
                    fig_name = sprintf('%s %s', str_title, strrep(coo.name, '_', ' '));
                    str_title = fig_name;
                else
                    name_list = [name_list {'UNKN'}];
                    fig_name = sprintf('Sigmas');
                end
                if ~isempty(coo.info.n_obs)
                    fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                    fh.UserData = struct('fig_name', fig_name);

                    if flag_prepro
                        ax(1) = subplot(2,1,1);
                        plotSep(coo.time.getMatlabTime, coo.info.s0_ip, '.-', 'MarkerSize', 10, 'LineWidth', 2); hold on;
                        ylabel('Sigma 0 PREPRO');
                        ax(2) = subplot(2,1,2);
                    end
                    plotSep(coo.time.getMatlabTime, coo.info.s0, '.-', 'MarkerSize', 10, 'LineWidth', 2); hold on;
                    ylabel('Sigma 0');
                    
                    yl(2) = max(yl(2), 2*perc(noNaN(coo.info.s0), 0.97));
                    xl(1) = min(xl(1), coo.time.first.getMatlabTime);
                    xl(2) = max(xl(2), coo.time.last.getMatlabTime);
                end
            end
            legend(name_list);
            if flag_prepro
                setAxis(fh,1);
                xlim(xl);
                setTimeTicks();
                title(sprintf('%s\\fontsize{5} \n', str_title));
                setAxis(fh,2);
                linkaxes(ax, 'x');
            end
            xlim(xl);
            ylim(yl);
            setTimeTicks();
            title(sprintf('%s\\fontsize{5} \n', str_title));
            
            fh_list = [fh_list fh];
            Core_UI.beautifyFig(fh);
            Core_UI.addBeautifyMenu(fh);
            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            
        end

        function fh_list = showSigmas_mr(coo_list, n_hours)
            % Simple display of sigma0 of the solutions
            %
            % SYNTAX
            %   fh_list = coo_list.showSigmas_mr();
            if nargin == 1
                n_hours = 24*4;
            end
                        
            [sigmas, time_sync] = coo_list.getProcessingSigmas(n_hours, true);
            [err_status, ~] = coo_list.getProcessingStatus(n_hours, true);
            
            sigmas = nan2zero(sigmas .* single(err_status == 0));
            time_sync = time_sync.getMatlabTime();
            
            fh = figure('Visible', 'off');
            fh.Name = sprintf('%03d: CooS0MR', fh.Number); fh.NumberTitle = 'off';
            fh_list = fh;
            fig_name = 'Coo_SigmasMR';
            fh.UserData = struct('fig_name', fig_name);
            
            for r = 1 : size(sigmas, 2)
                scatter(time_sync, r * ones(size(sigmas(:,r))), 30, double(sigmas(:,r))', 'filled'); hold on;
            end
            colormap(flipud(hot));

            yticks(1 : numel(coo_list));
            yticklabels(coo_list.getName());
            xlim([time_sync(1) time_sync(end)] + median(diff(time_sync), 'omitnan') .* [-1 1]);
            ylim([0.7 numel(coo_list) + 0.3]);
            setTimeTicks(28);
            grid on;
            title(sprintf('Epochs sigma0\\fontsize{5} \n'), 'FontName', 'Open Sans');
            Core_UI.beautifyFig(fh);
            
            Core_UI.addExportMenu(fh);
            Core_UI.addBeautifyMenu(fh);            
            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end
        
        function fh_list = showFixingRatio(coo_list)
            % Simple display of the fixing ratio
            %
            % SYNTAX
            %   fh_list = coo_list.showFixingRatio(n_hours);
            fh_list = [];
            coo_list = coo_list(not(coo_list.isEmpty));
            
            yl = [0 100];
            xl = [inf -inf];
            
            str_title = 'Fixing Ratio';
            if not(isempty(coo_list))
                fh = figure('Visible', 'off');  Core_UI.beautifyFig(fh);
            end
            name_list = {};
            for coo = coo_list(:)'
                if not(isempty(coo.name))
                    name_list = [name_list {coo.name}];
                    fig_name = sprintf('%s %s', str_title, coo.name);
                    str_title = fig_name;
                else
                    name_list = [name_list {'UNKN'}];
                    fig_name = sprintf('Fixing Ratio');
                end
                if ~isempty(coo.info.n_obs)
                    fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                    fh.UserData = struct('fig_name', fig_name);

                    plotSep(coo.time.getMatlabTime, coo.info.fixing_ratio, '.-', 'MarkerSize', 10, 'LineWidth', 2); hold on;
                    ylabel('Fixing ratio');
                    
                    xl(1) = min(xl(1), coo.time.first.getMatlabTime);
                    xl(2) = max(xl(2), coo.time.last.getMatlabTime);
                end
            end
            legend(name_list);
            xlim(xl);
            ylim(yl);
            setTimeTicks();
            title(sprintf('%s\\fontsize{5} \n', str_title));
            
            fh_list = [fh_list fh];
            Core_UI.beautifyFig(fh);
            Core_UI.addBeautifyMenu(fh);
            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end        
        
        function fh_list = showFixingRatio_mr(coo_list, n_hours)
            % Simple display of fixing ratio of the solutions
            %
            % SYNTAX
            %   fh_list = coo_list.showFixingRatio_mr();
            if nargin == 1
                n_hours = 24*4;
            end
                        
            [fixing_ratio, time_sync] = coo_list.getProcessingFixingRatio(n_hours, true);
            [err_status, ~] = coo_list.getProcessingStatus(n_hours, true);
            
            fixing_ratio = nan2zero(fixing_ratio .* single(err_status == 0));
            time_sync = time_sync.getMatlabTime();
            
            fh = figure('Visible', 'off');
            fh.Name = sprintf('%03d: CooFRMR', fh.Number); fh.NumberTitle = 'off';
            fh_list = fh;
            fig_name = 'Coo_FixRatioMR';
            fh.UserData = struct('fig_name', fig_name);
            
            for r = 1 : size(fixing_ratio, 2)
                scatter(time_sync, r * ones(size(fixing_ratio(:,r))), 30, double(fixing_ratio(:,r))', 'filled'); hold on;
            end
            colormap(flipud(Cmap.get('viridis')));

            yticks(1 : numel(coo_list));
            yticklabels(coo_list.getName());
            xlim([time_sync(1) time_sync(end)] + median(diff(time_sync), 'omitnan') .* [-1 1]);
            ylim([0.7 numel(coo_list) + 0.3]);
            setTimeTicks(28);
            grid on;
            colorbar;
            title(sprintf('Epoch fixing ratio\\fontsize{5} \n'), 'FontName', 'Open Sans');
            Core_UI.beautifyFig(fh);
            
            Core_UI.addExportMenu(fh);
            Core_UI.addBeautifyMenu(fh);            
            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
        end

        function fh_list = showErrCodes(coo_list, n_hours, flag_now)
            % Display error codes for data processing sessions over time.
            %
            % This function visualizes the error status of data processing sessions
            % using color-coded markers. It splits the plotting into chunks based on
            % a maximum number of stations allowed per figure.
            %
            % INPUT:
            %   coo_list  - List of sessions or data instances for which error codes need to be visualized.
            %   n_hours   - (Optional) Number of hours for which data needs to be considered. Default is 96 hours.
            %   flag_now  - (Optional) Flag to indicate if data should be processed until current time. Default is true.
            %
            % OUTPUT:
            %   fh_list   - Array of figure handles for each created figure.
            %
            % SYNTAX:
            %   fh_list = showErrCodes(coo_list, n_hours, flag_now)
            %
            % EXAMPLE:
            %   fh_list = showErrCodes(coo_list, 72, true);
            %

            max_n_sta = 50;
            % Simple display of the error codes per session
            if nargin == 1
                n_hours = 24*4;
            end

            if nargin == 3
                [err_status, time_sync] = coo_list.getProcessingStatus(n_hours, flag_now);
            else
                [err_status, time_sync] = coo_list.getProcessingStatus(n_hours, true);
            end

            good_perc = sum(err_status == 0)./size(err_status,1)*100;
            for i = 1 : numel(good_perc)
                fprintf('%s %.2f\n', coo_list(i).getNameV3, good_perc(i));
            end
            fprintf('--------------------------------\n %s %.2f\n', 'TOTAL    ', sum(err_status(:) == 0)./numel(err_status)*100);
            time_sync = time_sync.getDateTime();

            num_sta = numel(coo_list);
            num_figures = ceil(num_sta / max_n_sta);
            fh_list = gobjects(num_figures, 1); % Preallocate figure handles array

            for fig_idx = 1:num_figures
                start_idx = (fig_idx-1)*max_n_sta + 1;
                end_idx = min(fig_idx*max_n_sta, num_sta);

                fh = figure('Visible', 'off');
                fh.Name = sprintf('%03d: CooErr', fh.Number);
                fh.NumberTitle = 'off';
                fig_name = 'Coo_ErrCodes';
                fh.UserData = struct('fig_name', fig_name);

                for r = start_idx:end_idx
                    % All ok
                    plot(time_sync(err_status(:,r) == 0), (r-start_idx+1) * ones(sum(err_status(:,r) == 0), 1), 'ok', 'MarkerFaceColor', Core_UI.GREEN, 'MarkerEdgeColor', Core_UI.GREEN, 'MarkerSize', 6); hold on;
                    % Err code happened
                    plot(time_sync(err_status(:,r) > 0), (r-start_idx+1) * ones(sum(err_status(:,r) > 0), 1), 'ok', 'MarkerFaceColor', Core_UI.GREY, 'MarkerEdgeColor', Core_UI.GREY, 'MarkerSize', 6); hold on;
                    % only pre-processing
                    id_nopp = rem(err_status(:,r), 100) >= 10;
                    plot(time_sync(id_nopp), (r-start_idx+1) * ones(sum(id_nopp), 1), 'ok', 'MarkerFaceColor', Core_UI.ORANGE, 'MarkerEdgeColor', Core_UI.ORANGE, 'MarkerSize', 6); hold on;
                    % not even pre-processed
                    plot(time_sync(rem(err_status(:,r),10) == 1), (r-start_idx+1) * ones(sum(rem(err_status(:,r),10) == 1), 1), 'ok', 'MarkerFaceColor', Core_UI.RED, 'MarkerEdgeColor', Core_UI.RED, 'MarkerSize', 6); hold on;
                    % no PPP or NET
                    plot(time_sync(err_status(:,r) >= 100), (r-start_idx+1) * ones(sum(err_status(:,r) >= 100), 1), 'ok', 'MarkerEdgeColor', Core_UI.BLACK, 'MarkerFaceColor', 'none', 'LineWidth', 1, 'MarkerSize', 6); hold on;
                end

                yticks(1 : end_idx - start_idx + 1);
                yticklabels(coo_list(start_idx:end_idx).getName());

                xlim([time_sync(1) time_sync(end)] + median(diff(time_sync), 'omitnan') .* [-1 1]);
                ylim([0.7 (end_idx-start_idx+1) + 0.3]);

                % Dummy plots for legend
                h1 = plot(NaN,NaN,'ok', 'MarkerFaceColor', Core_UI.GREEN, 'MarkerEdgeColor', Core_UI.GREEN, 'MarkerSize', 6); hold on;
                h2 = plot(NaN,NaN,'ok', 'MarkerFaceColor', Core_UI.GREY, 'MarkerEdgeColor', Core_UI.GREY, 'MarkerSize', 6); hold on;
                h3 = plot(NaN,NaN,'ok', 'MarkerFaceColor', Core_UI.ORANGE, 'MarkerEdgeColor', Core_UI.ORANGE, 'MarkerSize', 6); hold on;
                h4 = plot(NaN,NaN,'ok', 'MarkerFaceColor', Core_UI.RED, 'MarkerEdgeColor', Core_UI.RED, 'MarkerSize', 6); hold on;
                h5 = plot(NaN,NaN,'ok', 'MarkerEdgeColor', Core_UI.BLACK, 'MarkerFaceColor', 'none', 'LineWidth', 1, 'MarkerSize', 6); hold on;

                % Add legend
                legend([h1, h2, h3, h4, h5], {'No errors', 'Some errors', 'Only pre-processed', 'Not even pre-processed', 'No PPP or NET'}, 'Location', 'best');

                grid on;
                title(sprintf('Bad epochs\\fontsize{5} \n'), 'FontName', 'Open Sans');
                Core_UI.beautifyFig(fh);

                Core_UI.addExportMenu(fh);
                Core_UI.addBeautifyMenu(fh);
                fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;

                fh_list(fig_idx) = fh;
            end

        end

        function [pos_time, t_max, pos_diff, pos_diff_model, pos_std, flag_time, pos_flags, time_ko] = getPosDiff(mode, coo, baricenter, coo_ref)
            % Extract the baseline or the absolute position of a coordinate object
            %
            % INPUT
            %   mode        'XYZ' , 'ENU'
            %   coo         coordinate object            
            %   baricenter  coordinate object of the observations (median?)
            %   coo_ref     coordinate object list of the reference/s
            %
            % SYNTAX
            %   [pos_time, t_max, pos_diff, pos_diff_model, pos_std, flag_time, pos_flags] = getPosDiff(mode, coo, baricenter, coo_ref)
            pos_time = [];
            pos_flags = [];
            log = Core.getLogger;
            flag_rem_baricenter = ~isempty(baricenter);
            if isempty(coo_ref)
                rate = coo.getRate;
                
                if strcmpi(mode, 'XYZ')
                    if flag_rem_baricenter
                        pos_diff = (coo.getXYZ - baricenter.getXYZ) * 1e3;
                    else
                        pos_diff = (coo.getXYZ - coo.getMedianPos.getXYZ) * 1e3;
                    end
                    if any(coo.xyz_model(:))
                        if flag_rem_baricenter
                            pos_diff_model = (coo.getXYZ('model') - baricenter.getXYZ) * 1e3;
                        else
                            pos_diff_model = (coo.getXYZ('model') - coo.getMedianPos.getXYZ) * 1e3;
                        end
                    else
                        pos_diff_model = [];
                    end
                    pos_std = coo.getStdXYZ .* 1e3;
                elseif strcmpi(mode, 'ENU')
                    if any(coo.xyz_model(:))
                        if flag_rem_baricenter
                            pos_diff_model = coo.getLocal(baricenter, 'model') * 1e3;
                        else
                            pos_diff_model = coo.getLocal(coo.getMedianPos, 'model') * 1e3;
                        end
                    else
                        pos_diff_model = [];
                    end
                    if flag_rem_baricenter
                        pos_diff = coo.getLocal(baricenter) * 1e3;
                    else
                        pos_diff = coo.getLocal(coo.getMedianPos) * 1e3;
                    end
                    pos_std = coo.getStdENU .* 1e3;
                end
                flag_time = true;
                if isa(coo.time, 'GPS_Time') && ~coo.time.isEmpty
                    pos_time = coo.time.getMatlabTime;
                    if numel(pos_time) < size(pos_diff,1)
                        log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more coordinates than times\n plotting only the positions with time'))
                        if numel(pos_diff_model) == numel(pos_diff)
                            pos_diff_model = pos_diff_model(1:numel(pos_time),:);
                        end
                        pos_diff = pos_diff(1:numel(pos_time),:);
                    elseif numel(pos_time) > size(pos_diff,1)
                        log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more times than coordinates\n plotting only the first positions'))
                        pos_time = pos_time(1:size(pos_diff,1),:);
                    end
                else
                    flag_time = false;
                    pos_time = (1 : size(pos_diff, 1))';
                end
                pos_flags = coo.info.flag;
                t_max = pos_time(end);
            else % plot baseline                
                pos_diff = [];
                if isa(coo.time, 'GPS_Time') && ~coo.time.isEmpty
                    new_ref_time = coo.time.getEpoch(~isnan(coo.xyz(:,1))).getCopy; % reference could have a different rate
                    [sync_time, id_sync] = getSyncedTime([new_ref_time; coo.time]);
                    % keep only epochs that are in common
                    id_ok = sum(~isnan(id_sync),2) == 2;
                    sync_time.remEpoch(~id_ok);
                    id_sync = id_sync(id_ok,:);
                    idx1 = id_sync(:,1);
                    idx2 = id_sync(:,2);
                    pos_time = sync_time.getMatlabTime;
                    if ~isempty(coo.info.flag)
                        pos_flags = coo.info.flag(idx2);
                    end
                    if strcmpi(mode, 'XYZ')
                        p_xyz = coo.getXYZ;
                        % Inter (nearest neighbour) the reference to the point of the relative
                        r_xyz_original = coo_ref.getXYZ;
                        r_xyz = nan(new_ref_time.length, 3); 
                        for c = 1:3
                            r_xyz(:,c) = interp1(coo_ref.time.getMatlabTime, r_xyz_original(:,c), new_ref_time.getMatlabTime, 'nearest', 'extrap');
                        end
                        pos_diff = (p_xyz(idx2,:) - r_xyz(idx1,:))*1e3;
                        pos_diff = bsxfun(@minus, pos_diff, median(pos_diff,1, 'omitnan'));
                        p_xyz_model = coo.getXYZ('model');
                        if numel(p_xyz_model) == numel(p_xyz) && any(p_xyz_model(:))
                            r_xyz_model_original = coo_ref.getXYZ(model);
                            r_xyz_model = nan(new_ref_time.length, 3); 
                            for c = 1:3
                                r_xyz_model(:,c) = interp1(coo_ref.time.getMatlabTime, r_xyz_model_original(:,c), new_ref_time.getMatlabTime, 'nearest', 'extrap');
                            end
                            pos_diff_model = (p_xyz_model(idx2,:) - r_xyz_model(idx1,:))*1e3;
                            pos_diff_model = bsxfun(@minus, pos_diff_model, median(pos_diff_model,1, 'omitnan'));
                        else
                            pos_diff_model = [];
                        end
                        % Compute formal std of the baseline
                        std1_original = coo_ref.getStdXYZ.^2;
                        std1 = nan(new_ref_time.length, 3); 
                        for c = 1:3
                            std1(:,c) = interp1(coo_ref.time.getMatlabTime, std1_original(:,c), new_ref_time.getMatlabTime, 'nearest', 'extrap');
                        end                        
                        std2 = coo.getStdXYZ.^2;
                        % covariance propagation: cov1 + cov2 + 2 * cross
                        % cross covariance is missing, propagation is incomplete
                        pos_std = sqrt(std1(idx1,:) + std2(idx2,:)) .* 1e3;
                    elseif strcmpi(mode, 'ENU')
                        p_xyz = coo.getXYZ;
                        r_xyz_original = coo_ref.getXYZ;
                        r_xyz = nan(new_ref_time.length, 3); 
                        for c = 1:3
                            if coo_ref.time.length <= 1
                                r_xyz(:,c) = r_xyz_original(:,c);
                            else
                                r_xyz(:,c) = interp1(coo_ref.time.getMatlabTime, r_xyz_original(:,c), new_ref_time.getMatlabTime, 'nearest', 'extrap');
                            end
                        end                        
                        pos_diff = Coordinates.cart2local(median(r_xyz,1,'omitnan'),p_xyz(idx2,:) - r_xyz(idx1,:) )*1e3;

                        p_xyz_model = coo.getXYZ('model');
                        if numel(p_xyz_model) == numel(p_xyz) && any(p_xyz_model(:))
                            r_xyz_model_original = coo_ref.getXYZ('model');
                            r_xyz_model = nan(new_ref_time.length, 3); 
                            for c = 1:3
                                r_xyz_model(:,c) = interp1(coo_ref.time.getMatlabTime, r_xyz_model_original(:,c), new_ref_time.getMatlabTime, 'nearest', 'extrap');
                            end
                            pos_diff_model = Coordinates.cart2local(median(r_xyz,1,'omitnan'),p_xyz_model(idx2,:) - r_xyz_model(idx1,:) )*1e3;
                        else
                            pos_diff_model = [];
                        end

                        % Compute formal std of the baseline
                        std1_original = coo_ref.getStdENU.^2;
                        std1 = nan(new_ref_time.length, 3); 
                        for c = 1:3
                            if coo_ref.time.length <= 1
                                std1(:,c) = std1_original(:,c);
                            else
                                std1(:,c) = interp1(coo_ref.time.getMatlabTime, std1_original(:,c), new_ref_time.getMatlabTime, 'nearest', 'extrap');
                            end
                            
                        end  
                        std2 = coo.getStdENU.^2;
                        % covariance propagation: cov1 + cov2 + 2 * cross
                        % cross covariance is missing, propagation is incomplete
                        pos_std = sqrt(std1(idx1,:) + std2(idx2,:)) .* 1e3;
                    end
                    t_max = new_ref_time.last.getMatlabTime;
                    flag_time = true;
                else
                    p_xyz = coo.getXYZ;
                    r_xyz = coo_ref.getXYZ;
                    if numel(r_xyz) == numel(p_xyz)
                        pos_diff = Coordinates.cart2local(median(r_xyz,1,'omitnan'),p_xyz - r_xyz)*1e3;
                        pos_std = [];
                        pos_time = pos_time(1:size(pos_diff,1),:);
                        flag_time = false;
                        pos_flags = coo.info.flag;                    
                    else
                        log.addError(sprintf('No time in coordinates and number off coordinates in ref different from coordinate in the second receiver'))
                    end
                    p_xyz_model = coo.getXYZ('model');
                    if numel(p_xyz_model) == numel(p_xyz) && numel(r_xyz_model) == numel(p_xyz) && any(p_xyz_model(:))
                        r_xyz_model = coo_ref.getXYZ(model);
                        pos_diff_model = Coordinates.cart2local(median(r_xyz,1,'omitnan'),p_xyz_model - r_xyz_model)*1e3;
                    else
                        pos_diff_model = [];
                    end
                    t_max = pos_time(end);
                end
                if any(pos_diff_model)
                    pos_diff_model = bsxfun(@minus, pos_diff_model, median(pos_diff,1,'omitnan'));
                end
                pos_diff = bsxfun(@minus, pos_diff, robAdj(pos_diff')');
            end

            if any(pos_std)
                pos_std(pos_std == 0) = 100e3; % no std => std set to 100m
            end

            % Find non linearly sampled coordinates (anomaly)
            rate = median(diff(pos_time(:)));
            time_ko = abs((mod(pos_time, rate) -  robAdj(mod(pos_time, rate)'))) > rate/100;
            if any(time_ko)
                Core.getLogger.addWarning(sprintf('Coordinate %s object contains %d/%d non regularly sampled values', coo.getNameV3, sum(time_ko), numel(time_ko)));
                pos_flags(time_ko) = [];
                pos_time(time_ko) = [];
                pos_diff(time_ko,:) = [];
                pos_std(time_ko,:) = [];
                if ~isempty(pos_diff_model)
                    pos_diff_model(time_ko,:) = [];
                end
            end
        end

        function fh_list = showCoordinates(mode, coo_list, coo_ref, n_obs, rot_angle_cw)
            % Plot ENU or XYZ coordinates
            %
            % SYNTAX
            %   this.showCoordinates(coo_list);
            
            if strcmpi(mode, 'XYZ')
                mode = 'XYZ';
                axis_label = {'ECEF X', 'ECEF Y', 'ECEF Z'};
            else
                mode = 'ENU';
                axis_label = {'East', 'North', 'Up'};
            end

            if nargin < 3
                coo_ref = [];
            end
            if nargin >= 4 && n_obs > 0
                time_start = coo_ref.getTime.getEpoch(max(1, coo_ref.length() - n_obs));
                time_start.addIntSeconds(-coo_ref.getTime.getRate());                
                for c = 1:numel(coo_list)
                    coo_list(c) = coo_list(c).getCopy;
                    coo_list(c).rem(GPS_Time('1970-01-07 00:00:00'), time_start);
                end
            end
            % Remove baricentric position instead of the median value;
            flag_rem_baricenter = false;
            flag_no_interp = true; % do not show interpolated value
            flag_add_status = false; % visualize pre-processed and un-preprocessed data with o or x
            flag_confidence = true;

            fh_list = [];
            thr = iif(numel(coo_list) > 16, 1, 0.8);
            flag_distr = false;
            flag_jump = false;
            
            str_title{2} = sprintf('STD (vs smoothed signal)');
            str_title{3} = sprintf('STD (vs smoothed signal)');
            log = Core.getLogger();
            fh = [];
            if flag_rem_baricenter
                xyz = nan(numel(coo_list),3);
                for i = 1 : numel(coo_list)
                    xyz(i,:) = coo_list(i).getMedianPos().getXYZ;
                end
                baricenter = Coordinates.fromXYZ(median(xyz,1,'omitnan')); 
            else
                baricenter = [];
            end
            ns = 0;   
            yl = ({[-0 0], [-0 0], [-0 0]});
            xl = {[+inf 0], [+inf 0], [+inf 0]};
            legend_str = {};
            for i = 1 : numel(coo_list)
                coo = coo_list(i);

                if ~coo.isEmpty && any(coo.xyz(:))

                    if ~isempty(coo_ref)
                        coo_ref = coo_ref.getMergedPos();
                    end
                    [pos_time, t_max, pos_diff, pos_diff_model, pos_std, flag_time, pos_flags, time_ko] = getPosDiff(mode, coo, baricenter, coo_ref);
                    time_id = find(~time_ko); time_id = time_id(logical(pos_flags));
                    if nargin >=5 && ~isempty(rot_angle_cw)
                        pos_diff = Coordinates.rotateDisplacements(pos_diff, -rot_angle_cw);
                        pos_diff_model = Coordinates.rotateDisplacements(pos_diff_model, -rot_angle_cw);
                        % ToDo pos_std = Coordinates.rotateDisplacementsSigmas(pos_std, rot_angle_cw);
                    end
                    if isempty(coo_ref)
                        rate = coo.getRate;
                        % single coordinates (no baseline)
                        if not(isempty(str_title{1}))
                            old_names = regexp(str_title{1}, '.*Position', 'match', 'once');
                            old_names = old_names(1:(end-9));   
                            old_val = regexp(str_title{1}, '(?<=signal\)).*', 'match', 'once');
                        else
                            old_names = '';
                            old_val = '';
                        end
                        if not(isempty(coo.name))
                            %str_title{1} = sprintf('%s %s\nPosition stability %s [mm] @ %gs\nSTD (vs smoothed signal)%s', old_names, coo.name, mode, rate, old_val);
                            str_title{1} = sprintf('Position stability %s [mm] @ %gs\nSTD (vs smoothed signal)%s', mode, rate, old_val);
                            fig_name = sprintf('d%s %s', mode, coo.name);
                            legend_str = [legend_str {coo.name}];
                        else
                            %str_title{1} = sprintf('%s Position stability %s [mm] @ %gs\nSTD (vs smoothed signal)%s', old_names, mode, rate, old_val);
                            str_title{1} = sprintf('Position stability %s [mm] @ %gs\nSTD (vs smoothed signal)%s', mode, rate, old_val);
                            fig_name = sprintf('d%s', mode);
                            legend_str = [legend_str {' - '}];
                        end
                    else % baseline mode
                        if not(isempty(coo.name))
                            rate = coo_ref.getRate;
                            if sum(~coo_list.isEmpty) > 1
                                ns = ns + 1;
                                if isempty(fh)
                                    %str_title{1} = sprintf('%s vs %s @ %gs\nBaseline stability %s \nSTD (vs smoothed signal)', coo.getName, coo_ref.name, rate, mode);
                                    str_title{1} = sprintf('Baseline stability %s vs %s [mm] @ %gs\nSTD (vs smoothed signal)', mode, coo_ref.name, rate);
                                    fig_name = sprintf('d%s %s vs %s', mode, coo.name, coo_ref.name);
                                    legend_str = [legend_str {coo.name}];
                                else
                                    offset = 7*(ns-1) - 3;
                                    % str_title{1} = [str_title{1}(1:offset) ' - ' coo.getName str_title{1}((offset+1):end)];
                                    fig_name = [fig_name(1:(offset + 2 + numel(mode))) ' - ' coo.getName fig_name((offset + 3 + numel(mode)):end)];
                                    legend_str = [legend_str {coo.name}];
                                end
                            else
                                [bsl, dh] = coo.getBaseline(coo_ref);
                                if bsl > 1e3
                                    str_title{1} = sprintf('%s - %s @ %gs\nBaseline stability %s (length %.2f Km, up %.2f m)\nSTD (vs smoothed signal)', coo.name, coo_ref.name, rate, mode, bsl*1e-3, dh);
                                else
                                    str_title{1} = sprintf('%s - %s @ %gs\nBaseline stability %s (length %.2f m, up %.2f m)\nSTD (vs smoothed signal)', coo.name, coo_ref.name, rate, mode, bsl, dh);
                                end
                                fig_name = sprintf('d%s %s - %s', mode, coo.name, coo_ref.name);
                                legend_str = [legend_str {coo.name}];

                            end
                        else
                            str_title{1} = sprintf('Baseline stability %s [mm]\nSTD (vs smoothed signal)', mode);
                            fig_name = sprintf('d%s bsl', mode);
                            legend_str = [legend_str {'UNKN'}];
                        end
                    end
                    if numel(coo_list) == 1 || isempty(fh)
                        fh = figure('Visible', 'off');
                        Core_UI.beautifyFig(fh);
                    end
                    if flag_distr
                        subplot(3, 12, 1:9);
                    else
                        subplot(3,1,1);
                    end
                    fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                    fh.UserData = struct('fig_name', fig_name);

                    if size(pos_diff, 1) > 1
                        if numel(coo_list) == 1
                            color_order = Core_UI.getColor(1:3,3);
                        else
                            color_order = Core_UI.getColor(i * [1 1 1], numel(coo_list));
                        end
                        
                        data = {};
                        data_component = {};
                        clear lh;
                        for c = 1 : 3
                            % MAIN PLOT --------------------------------------------------------------------------------------------------------
                            setAxis(fh, c);
                            if flag_distr
                                subplot(3, 12, (c-1)*12 + (1:9));
                            else
                                subplot(3, 1, c);
                            end
                            data_component{c} = pos_diff(:,c);
                            flag_detrend = false;
                            if flag_detrend
                                data_component{c}(abs(data_component{c}) > 100) = nan;
                                data_component{c} = detrend(data_component{c}, 'omitnan');
                            end
                            setAxis(fh, c);
                            
                            % Plot confidence level
                            if flag_confidence
                                if any(pos_std(:))
                                    yyaxis right; ylabel('Formal std [mm]');
                                    tmp_ax = gca;
                                    tmp_ax.YColor = min(1, color_order(c,:)+0.2);
                                    Core_Utils.patchSep(pos_time(:), pos_std(:, c), color_order(c,:), 'FaceColor', color_order(c,:),'EdgeColor', tmp_ax.YColor,'FaceAlpha',0.15, 'EdgeAlpha', 0.25, 'HandleVisibility','off'); hold on;
                                    ylim([0 min(1e2,max(0.5, 4 * perc(pos_std(:),0.8)))]);
                                    yyaxis left;
                                end
                            end
                            
                            if thr < 1
                                % If there was no PPP nor NET, consider the data as bad
                                lid_pro_ok = logical(coo(1).isPPPOk(pos_flags)) | logical(coo(1).isNETOk(pos_flags)) | logical(pos_flags == 0); % this last condition is needed for files imported from text with no flag information
                                if all(lid_pro_ok == 0)
                                    % If none is ok, maybe all are from juste pre-processing
                                    lid_pro_ok = coo(1).isPreProcessed(pos_flags);
                                else
                                    qi = coo(1).getQualityIndex > 0.2;
                                    lid_pro_ok = lid_pro_ok & qi(~time_ko);
                                end
                                if isempty(lid_pro_ok)
                                    lid_pro_ok = true(size(data_component{c}));
                                end
                                pos_var = nan2zero(pos_std(:, c).^2);
                                flag_var = any(pos_var);
                                pos_var(pos_var == 0) = 100.^2; % 100 meters of std

                                %[sync_time, id_sync] = getSyncedTime([coo.time; GPS_Time(pos_time)]);
                                %id_in_coo = ~isnan(id_sync(:,2));
                                %lid_pro_ok = lid_pro_ok(id_in_coo); % cut to the size of pos_data
                                not_nan =  ~isnan(data_component{c});
                                if numel(lid_pro_ok) == numel(not_nan)
                                    not_nan = not_nan & lid_pro_ok;
                                end
                                data_smooth =  nan(size(pos_time));
                                trend =  nan(size(pos_time));
                                lid_ko = true(size(pos_time));
                                data{c} = nan(size(data_component{c}));
                                if flag_var
                                    [data{c}(not_nan), lid_ko(not_nan), trend(not_nan), data_smooth(not_nan), jump_list] = Coordinates.cooFilter([pos_time(not_nan) data_component{c}(not_nan) pos_var(not_nan)], 0.8, 7);
                                else
                                    [data{c}(not_nan), lid_ko(not_nan), trend(not_nan), data_smooth(not_nan), jump_list] = Coordinates.cooFilter([pos_time(not_nan) data_component{c}(not_nan)], 0.9, 4);
                                end                                
                                %DEBUG: figure; plot(data_component{c}); hold on; plot(data{c},'o')

                                % Sum outliers to badly processed data
                                if numel(lid_pro_ok) > 1
                                    if any(lid_pro_ok)
                                        lid_ko = lid_ko | ~lid_pro_ok;
                                    end
                                    data{c}(~lid_pro_ok) = nan; % remove data with no pre-processing or ppp/net
                                end

                                % with more than 60% of outliers do not perform outlier detection
                                if sum(lid_ko(not_nan))/numel(lid_ko(not_nan)) > 0.6
                                    data{c} = data_component{c};
                                    lid_ko = ~lid_pro_ok;
                                end
                                if any(pos_diff_model(:))
                                    data_smooth = pos_diff_model(:,c);
                                end
                                setAxis(fh, c);
                                
                                % Plot all the data
                                Core_Utils.plotSep(pos_time, data_component{c}, '.', 'MarkerSize', 7, 'LineWidth', 1.6, 'Color', [0.5 0.5 0.5 0.4]);
                                if flag_add_status && ~isempty(pos_flags)
                                    % data with only pre-processing
                                    tmp = data_component{c}; tmp(~(coo(1).isPreProcessed(pos_flags) & ~(logical(coo(1).isPPPOk(pos_flags)) | logical(coo(1).isNETOk(pos_flags))))) = nan;
                                    if any(tmp)
                                        Core_Utils.plotSep(pos_time, tmp, 'o', 'MarkerSize', 5, 'Color', min(1, color_order(c,:) + 0.2));
                                    end
                                    % data not even pre-processed
                                    tmp = data_component{c}; tmp(logical(coo(1).isPreProcessed(pos_flags))) = nan;
                                    if any(tmp)
                                        Core_Utils.plotSep(pos_time, tmp, 'x', 'MarkerSize', 5, 'Color', min(1, color_order(c,:) + 0.2));
                                    end
                                end

                                % Plot only supposely good data
                                Core_Utils.plotSep(pos_time, data{c}, '.-', 'MarkerSize', 6, 'LineWidth', 1.6, 'Color', color_order(c,:));
                                
                                % Plot smoothed signal
                                if ~flag_no_interp
                                    if any(data_smooth) && (std(data_smooth, 'omitnan') < 2*std(data{c}, 'omitnan') || std(data_smooth, 'omitnan'))
                                        jump_list = [jump_list, sum(not_nan)];
                                        for j = 1 : numel(jump_list)-1
                                            plot(pos_time(jump_list(j)+1:jump_list(j+1)), data_smooth(jump_list(j)+1:jump_list(j+1)), '-.', 'MarkerSize', 5, 'LineWidth', 1, 'Color', [max(0, color_order(c,:)-0.3) 0.8]);
                                        end
                                    else
                                        Core.getLogger.addWarning('Coordinate model not displaied, it is too bad');
                                    end
                                end
                                
                                ylj = [-30 30]*1e3;
                                tj = pos_time(not_nan);
                                if flag_jump
                                    for j = 2 : numel(jump_list)
                                        plot(tj([jump_list(j) jump_list(j)]) + (rate/2)/86400 , ylj, '-.', 'LineWidth', 2, 'Color', [0 0 0 0.75]); hold on;
                                    end
                                end
                                
                            else
                                data{c} = data_component{c};
                                setAxis(fh, c);
                                % Plot data
                                Core_Utils.plotSep(pos_time, data_component{c}, '.-', 'MarkerSize', 8.5, 'LineWidth', 1.6, 'Color', color_order(c,:));
                                lid_ko = false(numel(pos_time), 1);
                                
                                if any(pos_diff_model(:))
                                    data_smooth = pos_diff_model(:,c);
                                else
                                    % Use a simple trend on data as smothed signal
                                    data_smooth = Core_Utils.interp1LS(pos_time(~isnan(pos_diff(:,1))), pos_diff(~isnan(pos_diff(:,c)),c), 1, pos_time);
                                end
                            end                            
                            setAxis(fh, c);
                            ax(4-c) = gca(fh);
                            if (pos_time(end) > pos_time(1))
                                xl_tmp = [pos_time(1) t_max];
                                xl{c} = [min([xl{c}(1), xl_tmp(1)]) max([xl{c}(2), xl_tmp(2)])];
                                xlim(xl{c});
                            end
                            yl_tmp = minMax(data{c}) + max(2, round(0.1 * diff(minMax(data{c})))) * [-1 1];
                            yl{c} = [min([-20, yl{c}(1), yl_tmp(1)]) max([20, yl{c}(2), yl_tmp(2)])];
                            ylim(yl{c});
                            h = ylabel([axis_label{c} ' [mm]']); h.FontWeight = 'bold';
                            grid on;
                            
                            %str_title{c} = sprintf('%s %s%.2f', str_title{c}, iif(i>1, '- ', ''), std((data{c}(~lid_ko) - data_smooth(~lid_ko)), 'omitnan'));
                            %h = title(str_title{c}, 'interpreter', 'none'); h.FontWeight = 'bold';
                            legend_str{end} = sprintf('%s | %.2f', legend_str{end}, std((data{c}(~lid_ko) - data_smooth(~lid_ko)), 'omitnan'));                            
                        end
                        % Only lines are legend elements
                        for c = 1 : 3
                            ax_tmp = setAxis(fh, c);
                            i = 0;
                            k = 0;
                            legend_cell = {};
                            for j = numel(ax_tmp.Children) : -1 : 1 % build the elegend elements
                                try
                                    if ax_tmp.Children(j).Type(1) == 'l' % its a line
                                        k = k + 1;
                                        if strcmp(ax_tmp.Children(j).LineStyle, '-')
                                            i = i + 1;
                                            legend_cell{k} =  legend_str{i};
                                        else
                                            legend_cell{k} = '';
                                        end
                                    end
                                catch
                                    legend_cell{k} = '';
                                end
                            end
                            lh(c) = legend(legend_cell, 'location', 'northeastoutside');
                        end
                        sgtitle(str_title{1});
                        set(lh(2:3), 'Visible','off');
                        
                        if flag_distr
                            for c = 1 : 3
                                % DISTRIBUTION PLOT ------------------------------------------------------------------------------------------------
                                
                                subplot(3, 12, c*12 + (-1:0));
                                %%
                                % Get distribution of data
                                
                                [n_res, x] = hist(data{c}, max(round(diff(minMax(data{c}))/1), 20));
                                rate_x = mean(diff(x));
                                padding = 1:0.5:10;
                                x = [-fliplr(padding) * rate_x + min(x),  x, padding * rate_x + max(x)];
                                n_res = [0*padding n_res 0*padding];
                                n_res_var = 1./ n_res.^2;
                                
                                % Normalize to one
                                n_res = n_res ./ sum(n_res);
                                
                                % Get distribution limits
                                lim = minMax(x);
                                x_out = linspace(lim(1), lim(end), 1000);
                                
                                [~, ~, ~, y] = splinerMat(x', n_res', max(3*rate_x, 2), 1e-4, x_out);
                                %x_out = x;
                                %y = n_res';
                                
                                ax_tmp = setAxis(fh, c + 3);
                                patch([max(0,y); 0; 0], [x_out, x_out(end), x_out(1)]', min(1, color_order(c,:) + 0.2), ...
                                    'FaceColor', min(1, color_order(c,:) + 0.2), ...
                                    'EdgeColor', color_order(c,:), ...
                                    'FaceAlpha', 0.2, ...
                                    'EdgeAlpha', 1,  ...
                                    'HandleVisibility', 'off'); 
                                ax_tmp.XAxis.TickValues = [];
                                ylim([min(x(1), yl{c}(1)) max(x(end), yl{c}(2))]);
                                xlim([-0.01 max(y)+0.1]);
                                
                                title('Distribution');
                            end
                        end
                        linkaxes(ax, 'x');
                        grid on;
                    else
                        log.addMessage('Plotting a single point static coordinates is not yet supported');
                    end
                    fh_list = unique([fh_list fh]);
                end
            end
            for f = 1:numel(fh_list)
                Core_UI.beautifyFig(fh_list(f));
                Core_UI.addBeautifyMenu(fh_list(f));
                Core_UI.addLineMenu(fh_list(f));

                fh_list(f).Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                if flag_time
                    for a = 1 : 3
                        setAxis(fh_list(f), a);
                        setTimeTicks(6);
                    end
                end
            end
        end
        
        function fh = showCoordinatesENU(coo_list, coo_ref, n_obs, rot_angle_cw)
            % Plot East North Up coordinates
            %
            % SYNTAX
            %   this.showCoordinatesENU(coo_list, coo_ref,n_obs);
            
            switch nargin
                case 1, fh = showCoordinates('ENU', coo_list);
                case 2, fh = showCoordinates('ENU', coo_list, coo_ref);
                case 3, fh = showCoordinates('ENU', coo_list, coo_ref, n_obs);
                case 4, fh = showCoordinates('ENU', coo_list, coo_ref, n_obs, rot_angle_cw);
            end
        end
        
        function fh = showPositionXYZ(coo_list)
            % Plot X Y Z coordinates
            %
            % SYNTAX
            %   this.showPositionXYZ(coo_list);
            fh = showCoordinatesXYZ(coo_list);
        end
        
        function fh = showCoordinatesXYZ(coo_list, coo_ref, n_obs)
            % Plot X Y Z coordinates
            %
            % SYNTAX
            %   this.showCoordinatesXYZ(coo_list, coo_ref,n_obs);
            switch nargin
                case 1, fh = showCoordinates('XYZ', coo_list);
                case 2, fh = showCoordinates('XYZ', coo_list, coo_ref);
                case 3, fh = showCoordinates('XYZ', coo_list, coo_ref, n_obs);
            end
        end
        
        function fh = showPositionPlanarUp(coo_list)
            % Plot East North Up coordinates
            %
            % SYNTAX
            %   this.showPositionENU(coo_list);
            fh = showCoordinatesPlanarUp(coo_list);
        end
        
        function fh_list = showCoordinatesPlanarUp(coo_list, id_ref, n_obs)
            % Plot East North Up coordinates
            %
            % SYNTAX
            %   this.showCoordinatesPlanarUp(coo_list);
            
            if (nargin < 2) || isempty(id_ref)
                id_ref = [];
            else
                if ischar(id_ref)
                    [~, id_ref] = coo_list.get(id_ref);
                elseif isa(id_ref, 'Coordinates')
                    coo_ref = id_ref;
                    [~, id_ref] = coo_list.get(coo_ref.getName);
                    if isempty(id_ref)
                        coo_list = [coo_list; coo_ref];
                        id_ref = numel(coo_list);
                    end
                end
            end

            coo_ref = coo_list(id_ref);           
            % Remove baricentric position instead of the median value;
            flag_rem_baricenter = false;
            flag_no_interp = true;

            if nargin > 2 && n_obs > 0
                time_start = coo_ref.getTime.getEpoch(max(1, coo_ref.length() - n_obs));
                time_start.addIntSeconds(-coo_ref.getTime.getRate());
                for c = 1:numel(coo_list)
                    coo_list(c) = coo_list(c).getCopy;
                    coo_list(c).rem(GPS_Time('1970-01-07 00:00:00'), time_start);
                end
            end

            fh_list = [];
            thr = iif(numel(coo_list) > 16, 1, 0.8);
            flag_distr = false;
            flag_jump = false;
            
            log = Core.getLogger();
            fh_list = []; fh = [];
            if flag_rem_baricenter
                xyz = nan(numel(coo_list),3);
                for i = 1 : numel(coo_list)
                    xyz(i,:) = coo_list(i).getMedianPos().getXYZ;
                end
                baricenter = Coordinates.fromXYZ(median(xyz,1,'omitnan')); 
            else
                baricenter = [];
            end

            flag_multirec = sum(~coo_list.isEmpty) > 1;
            str_legend_prob = {};
            str_legend = {};
            for i = 1 : numel(coo_list)
                coo = coo_list(i);

                if ~coo.isEmpty && any(coo.xyz(:))
                    
                    if ~isempty(coo_ref)
                        coo_ref = coo_ref.getMergedPos();
                    end

                    [pos_time, t_max, enu_diff, pos_diff_model, pos_std, flag_time] = getPosDiff('ENU', coo, baricenter, coo_ref);
                    if nargin == 1 || isempty(coo_ref)
                        if isa(coo.time, 'GPS_Time') && ~coo.time.isEmpty
                            if numel(pos_time) < size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more coordinates than times\n plotting only the positions with time'))
                            elseif numel(pos_time) > size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more times than coordinates\n plotting only the first positions'))
                            end
                        end
                    end
                    
                    if size(enu_diff, 1) > 1
                        if isempty(fh)
                            str_title{1} = sprintf('%sPosition detrended planar Up [mm]\nSTD (vs smoothed signal)', iif(flag_multirec, '', [coo.getName ' ']));
                            str_title{3} = sprintf('STD ENU (vs smoothed signal)');
                            fh = figure('Visible', 'off'); Core_UI.beautifyFig(fh);
                            fh_list = [fh_list; fh];
                        end
                        if flag_multirec
                            fh.Name = sprintf('%03d: dPUP MR', fh.Number); fh.NumberTitle = 'off';
                            color_order = Core_UI.getColor(i * [1 1 1], numel(coo_list));
                        else
                            fh.Name = sprintf('%03d: dPUP', fh.Number); fh.NumberTitle = 'off';
                            color_order = Core_UI.getColor(1:3,3);
                        end
                                                
                        % OUTLIER DETECTION ------------------------------
                        trend =  nan(size(pos_time,1),3);
                        for c = 1:3
                            pos_var = nan2zero(pos_std(:, c).^2);
                            flag_var = any(pos_var);
                            pos_var(pos_var == 0) = 100.^2; % 100 meters of std
                            not_nan =  ~isnan(enu_diff(:,c));
                            data_smooth =  nan(size(pos_time));                            
                            lid_ko = true(size(pos_time));
                            trend_c =  nan(size(pos_time));
                            data{c} = nan(size(enu_diff(:,c)));
                            if flag_var
                                [data{c}(not_nan), lid_ko(not_nan), trend_c(not_nan), data_smooth(not_nan), jump_list] = Coordinates.cooFilter([pos_time(not_nan) enu_diff(not_nan,c) pos_var(not_nan)], 0.8, 7);
                            else
                                [data{c}(not_nan), lid_ko(not_nan), trend_c(not_nan), data_smooth(not_nan), jump_list] = Coordinates.cooFilter([pos_time(not_nan) enu_diff(not_nan,c)], 0.8, 7);
                            end

                            trend(:,c) = Core_Utils.interp1LS(pos_time(~isnan(data{c})), enu_diff(~isnan(data{c}),c), 1, pos_time);
                            enu_diff(:,c) = enu_diff(:,c) - trend(:,c);

                            % with more than 60% of outliers do not perform outlier detection
                            if sum(lid_ko(not_nan))/numel(lid_ko(not_nan)) > 0.6
                                data{c} = enu_diff(:,c);
                                lid_ko = lid_ko & false;
                            end
                            if any(pos_diff_model(:))
                                data_smooth = pos_diff_model(:,c);
                            end
                        end
                        
                        Core_UI.beautifyFig(fh);
                        fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                        drawnow
                        
                        % Plot parallel
                        max_e = ceil(max(abs(minMax(data{1})))/5) * 5;
                        max_n = ceil(max(abs(minMax(data{2})))/5) * 5;
                        max_r = ceil(sqrt(max_e^2 + max_n^2) / 5) * 5;
                        if max_r == 0
                            if ~flag_multirec
                                delete(fh);
                                fh_list(end) = [];
                                Core.getLogger.addWarning(sprintf('Coordinates are fixed or not available for "%s"', coo.getNameV3));
                            end
                        else
                            %% Polar plot -------------------------------------
                            fh.Visible = 'off';
                            subplot(3,2,[1 3]);
                            
                            % Plot circles of precision
                            az_l = 0 : pi/200: 2*pi;
                            % dashed
                            id_dashed = serialize(bsxfun(@plus, repmat((0:20:395)',1,5), (1:5)));
                            az_l(id_dashed) = nan;
                            step = max(1, floor((max_r/10)/5)) * 10;
                            if step > 10
                                step = round(step/50)*50;
                            end
                            if step > 100
                                step = round(step/100)*100;
                            end
                            decl_s = ((step : step  : max_r));
                            for d = decl_s
                                x = cos(az_l).*d;
                                y = sin(az_l).*d;
                                plot(x,y,'color',[0.6 0.6 0.6], 'LineWidth', 2); hold on;
                                x = cos(az_l).*(d-step/2);
                                y = sin(az_l).*(d-step/2);
                                plot(x,y,'color',[0.75 0.75 0.75], 'LineWidth', 2); hold on;
                                str_legend = [str_legend {'', ''}];
                            end
                            
                            plot(data{1} + trend(:,1), data{2} + trend(:,2), 'o', 'MarkerSize', 4, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                            axis equal;
                            h = ylabel('East [mm]'); h.FontWeight = 'bold';
                            h = xlabel('North [mm]'); h.FontWeight = 'bold';
                            ylim(max_r * [-1 1]);
                            xlim(max_r * [-1 1]);
                            grid on;
                            
                            str_legend = [str_legend {sprintf('STD ENU %.2f - %.2f - %.2f', std(data{1}, 'omitnan'), std(data{2}, 'omitnan'), std(data{3}, 'omitnan'))}];
                            legend(str_legend);

                            h = title(str_title{1}, 'interpreter', 'none'); h.FontWeight = 'bold';
                            h.FontWeight = 'bold';
                            
                            %% Distribution plot ---------------------------
                            subplot(3,2,[2 4]);
                            planar_full = sqrt(enu_diff(:,1).^2 + enu_diff(:,2).^2);
                            planar_filt = sqrt(data{1}.^2 + data{2}.^2);
                            
                            sigma1 = std(planar_filt, 'omitnan');
                            classes = linspace(0,6*sigma1, 600);
                            if flag_multirec
                                color = color_order(1,:);
                            else
                                c = hist(planar_filt, classes);
                                c_sum = cumsum(c)./numel(planar_full);
                                plot(classes, c_sum, 'LineWidth', 2); hold on;
                                color = [0 0 0];
                            end

                            c = hist(planar_full, classes);
                            c_sum = cumsum(c)./numel(planar_full);
                            if flag_multirec
                                plot(classes, c_sum, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                            else
                                plot(classes, c_sum, 'LineWidth', 2); hold on;
                            end

                            plot(1*[0; sigma1; sigma1], [c_sum(100); c_sum(100); 0], '--', 'Color', [color 0.5]);
                            plot(2*[0; sigma1; sigma1], [c_sum(200); c_sum(200); 0], '--', 'Color', [color 0.5]);
                            plot(3*[0; sigma1; sigma1], [c_sum(300); c_sum(300); 0], '--', 'Color', [color 0.5]);

                            grid minor;
                            ylabel('Probability');
                            xlabel('Error [mm]');
                            if flag_multirec
                                title(sprintf('Probability distribution - planar', sigma1));
                                str_legend_prob = [str_legend_prob {sprintf('%s STD %.2f', coo.getName, sigma1), '', '', ''}];
                                legend(str_legend_prob, 'Location', 'SouthEast');
                            else
                                title(sprintf('Probability distribution - planar STD: %.2f', sigma1));
                                legend('without outliers', 'including outliers', 'Location', 'SouthEast');
                            end

                            %% UP plot -------------------------------------
                            
                            set(0, 'CurrentFigure', fh);
                            subplot(3,2,[5 6]);
                            up = data{3};
                            Core_Utils.plotSep(pos_time, up, '.-', 'MarkerSize', 8.5, 'LineWidth', 1.6, 'Color', color_order(3,:)); hold on;
                            ax(1) = gca(fh);
                            if (pos_time(end) > pos_time(1))
                                xlim([pos_time(1) pos_time(end)]);
                            end
                            yl = [-1 1] * max(abs([perc(up, 0.05)*3 perc(up, 0.95)*3]));
                            ylim([min(-20, yl(1)) max(20, yl(2))]);
                            setTimeTicks(); h = ylabel('Up [mm]'); h.FontWeight = 'bold';
                            grid on;
                            drawnow;
                        end
                    else
                        log.addMessage('Plotting a single point static coordinates is not yet supported');
                    end
                end
            end
            for fh = fh_list(:)'
                Core_UI.beautifyFig(fh);
                Core_UI.addBeautifyMenu(fh);
                fh.Visible = iif(Core_UI.isHideFig, 'off', 'on');
            end
        end
        
        function [fh, mask] = showSkyMask(coo, varargin)
            % Get a predicted sky mask from DTM
            % Warning it does not handle yet areas close to 180 -180 longitudes
            %
            % SYNTAX
            %   [mask, fh] = coo.showSkyMask(...)
            %
            % INPUT
            %   Input are defined as a list of (couples option, value)
            %   e.g. coo.getSkyMask('rec', rec, 'flag:_fig', true, ...)
            %
            %   'rec'          use an external receiver
            %                      GNSS_Station
            %                      default:  empty
            %
            %   'use_dtm_h0'   use the height of the DTM as height of the
            %                  point, if folse use the height in rec/coo
            %                      flag:     true / false
            %                      default:  true
            %
            %   'h_antenna'    additional heigth of the antenna
            %                      numeric
            %                      default:  0
            %
            %   'h_vegetation'  additional heigth of the DTM
            %                      numeric
            %                      default:  0
            %
            %   'az_mask'      array of azimuth for the mask
            %                      numeric array [deg]
            %                      default:  (0:0.5:360)'
            %
            %   'az_res'       when az_mask is of one element use this number as resolution5295
            %                      numeric [deg]
            %                      default:  0.5
            %
            %   'min_el'       minimum elevation to display [deg]
            %                      numeric [deg]
            %                      default:  -7.5
            %
            %  - FIGURE ---------------------------------------------------
            %
            %   'name'         display a different name of the station
            %                      char
            %                      default:  name of coordinates
            %
            %   'marker'       display a different marker of the station on the map
            %                      char
            %                      default:  marker of coordinates
            %
            %  - ORBITS ---------------------------------------------------
            %    Note: orbits are computed using final GFZ products 
            %    (or the loaded core if a rec is passed)
            %
            %   'date_start'   first epoch to compute orbits
            %                      GPS_Time
            %                      default:  GPS_Time('2021-11-1 00:00:00')
            %
            %   'date_stop'    last epoch to compute orbits
            %                      GPS_Time
            %                      default:  GPS_Time('2021-11-1 23:59:59')
            % 
            %   'sys_list'     array of constellation to use
            %                      logical [7x1] || char
            %                      default:  'GRE'
            %             
            varargin = [{'show_fig', 'true'}, varargin];
            coo = coo.getMedianPos;
            [mask, fh] = coo.getSkyMask(varargin{:});
        end
        
        function [fh] = showMap(coo, varargin)
            % Display a map of the coordinates
            % Warning it does not handle yet areas close to 180 -180 longitudes
            %
            % SYNTAX
            %   fh = coo.showMap(...)
            %
            % INPUT
            %   Input are defined as a list of (couples option, value)
            %   e.g. coo.showMap('bg_type', 'map', 'id_ref', 5, ...)
            %
            %   'new_fig'          open a new figure or uses
            %                      flag:     true / false
            %                      default:  true
            %
            %   'fig_handle'       use this fig handle instead of creating a new one
            %                      default none
            %
            %   'use_model'        use modelled coordinates instead of observed
            %                      flag:     true / false
            %                      default:  false
            %
            %  - POINTS ---------------------------------------------------
            %
            %   'N_MAX_POINTS'     max number labels to display
            %                      integer:  0..inf
            %                      default:  50
            %
            %   'flag_tooltip'     display tooltip
            %                      flag:     true / false
            %                      default:  true
            %
            %   'id_ref'           id of the reference station
            %                      integer:  empty or 1..n_coo
            %                      default:  empty (no reference)
            %
            %   'label_size'       size of the labels
            %                      float:  1..n
            %                      default:  12
            %
            %   'flag_label'       display labels
            %                      flag:     true / false
            %                      default:  true
            %
            %   'flag_label_bg'    display label background (box)
            %                      flag:     true / false
            %                      default:  true
            %                      values are forced for some bg_type (false | map) or (true | dtm)
            %
            %   'label_color'      display label background (box)
            %                      rgb:      [1 1 1]
            %                      default:  white
            %
            %   'point_size'       Dimension of the points
            %                      number:   [1 x 1]
            %                      default:  10
            %
            %   'point_color'      Color of the points [RGB]
            %                      flag:     [3 x 1]
            %                      default:  Core_UI.LBLUE
            %
            %  - MAP TYPE -------------------------------------------------
            %
            %   'proj_type'        define the type of projection
            %                      flag:     'none' / 'equidistant' / 'UTM'
            %                      default:  'Equidistant'
            %
            %   'bg_type'          define the type of map background to download
            %                      (google maps, shaded DTM)
            %                      flag:     'none' / 'map' / 'dtm'
            %                      default:  'map'
            %
            %   'map_type'         see addMap function for the possible values
            %                      default: 'ArcGIS'
            %
            %   'dtm_resolution'   resolution of the DTM to download
            %                      flag:     '1' / '3' / 'high' / 'low'
            %                      default:  '1'
            %
            %   'shape'            global shapefile to use
            %                      flag:     'none' / 'coast' / 'fill' / '10m' / '30m' / '50m'
            %                      default:  'none'
            %
            %   'add_margin'       add further border to the points [deg]
            %                      array:    [2 x 1]
            %                      default:  [0 0]
            %
            %   'contour_lines'    add contour on the map
            %                      if numeric specify the distance in [m] between curves
            %                      flag:     true / false
            %                      default:  false
            %
            %  - VELOCITY ARROWS ------------------------------------------
            %
            %   'arrow_type'       type of arrow array to display, if none no arrows will be displayed
            %                      flag:     'none' / 'planar' / 'up'
            %                      default:  'none'
            %
            %   'use_disp'         show displacements instead of velocities (cumulative displacement till day 0)
            %                      flag:     true / false
            %                      default:  false
            %
            %   'animate'          animate_arrows
            %                      flag:     true / false
            %                      default:  false
            %
            %   'ani_export'       export animation
            %                      flag:     'none' / 'gif' / 'video'
            %                      default:  'none'
            %
            %   'ani_file_name'    animation file name
            %                      path:     '' / <full_path>
            %                      default:  '' == it will be automatic defined
            %
            %   'recompute_vel'    flag to recompute velocity vector from linear regrassion of the stored positions
            %                      if false use the velocity already stored within the Coordinate object
            %                      flag:     true / false
            %                      default:  true
            %
            %   'vel_lim'          max and min size of the velocity arrows [mm],
            %                      the biggest arrow is displayed with 1/6 of the size of the map
            %                      matrix:   [2 x 2] [max_planar min_planar;
            %                                         max_up     min_up]
            %                      default:  [100 10;
            %                                 100 10]
            %
            %   'cvel_lim'         caxis min and max for the arrow colormap [mm]
            %                      array:    [2 x 1]
            %                      default:  [60 60]
            %
            %   'cmap'             RGB colormap for the arrows
            %                      array:    [n_colors x 3]
            %                      default:  Cmap.get('c51'); cmap = flipud(cmap(20:36,:));
            
            coo = coo(not(coo.isEmpty));
            if not(isempty(coo))
                % Default args --------------------------------------------
                
                new_fig = true;
                use_model = false;
                coo_type = 'obs'; % this is later commanded by use_model, if true this field will be 'model'
                proj_type = 'equidistant';
                bg_type = 'map';          % it can be: none, map, dtm
                map_type = 'ArcGIS';
                flag_no_prj = false;
                dtm_resolution = 'auto';   % it can be: "low", "high"
                shape = 'none';            % it can be: none / coast / fill / 10m / 30m / 50m
                contour_lines = false;
                add_margin = [0 0];
                
                id_ref = [];
                N_MAX_POINTS = 50; % Max number of stations to use larger points / labels
                
                flag_tooltip = true;
                flag_label = true;
                flag_label_bg = false;
                flag_fix_label = true;
                label_color = [0.1 0.1 0.1];
                label_size = 12;
                point_size = 10;
                point_color = Core_UI.LBLUE;
                
                arrow_type = 'none';       % It can be: none / planar / none
                use_displacements = false;
                animate = false;
                ani_export = 'none';
                ani_file_name = '';
                recompute_vel = true;
                vel_lim = [100 50; 100 50];% scaling for arrows [mm] [max_x, min_x; max_y, min_y]
                cvel_lim = [60 60];        % max red alert in mm / h [planar up]
                cmap = Cmap.get('c51'); cmap = flipud(cmap(20:36,:));

                fill_val = [];
                fill_lim = [];

                fig_handle = [];

                % Parse args ----------------------------------------------
                
                force_label = false;
                args = varargin;
                if iscell(args) && numel(args) > 0
                    if iscell(args{1})
                        args = args{1};
                    end
                end
                if nargin > 1
                    a = 1;
                    while a < numel(args)
                        if ischar(args{a})
                            a = a + 1;
                            switch args{a-1}
                                case 'fig_handle'
                                    fig_handle = args{a};
                                case 'fill_val'
                                    fill_val = args{a};
                                case 'fill_lim'
                                    fill_lim = args{a};
                                case 'new_fig'
                                    new_fig = logical(args{a});
                                case 'use_model'
                                    use_model = logical(args{a});
                                case {'proj_type', 'type', 'proj'}
                                    if ismember(args{a}, {'UTM', 'equidistant', 'mercator'})
                                        proj_type = args{a};
                                    else
                                        proj_type = 'none';
                                    end
                                case 'bg_type'
                                    if ismember(args{a}, {'map', 'dtm'})
                                        bg_type = args{a};
                                    else
                                        bg_type = 'none';
                                    end
                                case {'map_type', 'provider'}
                                    map_type = args{a};
                                case 'dtm_resolution'
                                    if ismember(args{a}, {'low', 'high', '3', '1'})
                                        dtm_resolution = args{a};
                                    else
                                        dtm_resolution = 'auto';
                                    end
                                case 'shape'
                                    if ismember(args{a}, {'coast', 'fill', '10m', '30m', '50m'})
                                        shape = args{a};
                                    else
                                        shape = 'none';
                                    end
                                case 'arrow_type'
                                    if ismember(args{a}, {'planar', 'up'})
                                        arrow_type = args{a};
                                    else
                                        arrow_type = 'none';
                                    end
                                case 'use_disp'
                                    use_displacements = logical(args{a});
                                case 'animate'
                                    animate = logical(args{a});
                                case 'ani_export'
                                    ani_export = args{a};
                                case 'ani_file_name'
                                    ani_file_name = args{a};
                                case 'flag_recompute_vel'
                                    recompute_vel = logical(args{a});
                                case 'vel_lim'
                                    vel_lim = args{a};
                                case 'cvel_lim'
                                    cvel_lim = args{a};
                                case 'N_MAX_POINTS'
                                    N_MAX_POINTS = args{a};
                                case 'id_ref'
                                    id_ref = args{a};
                                case 'add_margin'
                                    add_margin = args{a};
                                case 'flag_tooltip'
                                    flag_tooltip = logical(args{a});
                                case 'flag_label'
                                    flag_label = logical(args{a}); force_label = true;
                                case 'label_size'
                                    label_size = args{a};
                                case 'label_color'
                                    label_color = args{a};
                                case 'flag_label_bg'
                                    flag_label_bg = args{a};
                                case 'flag_fix_label'
                                    flag_fix_label = args{a};
                                case 'contour_lines'
                                    contour_lines = args{a};
                                case 'point_size'
                                    point_size = args{a};
                                case 'point_color'
                                    point_color = args{a};
                                case 'cmap'
                                    cmap = args{a};
                                otherwise
                                    if ischar(args{a-1})
                                        log = Core.getLogger;
                                        log.addWarning(sprintf('Parameter "%s" unrecognized in Coordinates.showMap()', args{a-1}));
                                    end
                            end
                        end
                        a = a + 1;
                    end
                end
                
                id_ref(id_ref > numel(coo)) = [];
                
                if animate
                    animate =  animate && not(arrow_type(1) == 'n');
                    use_displacements = true;
                end
                
                if use_model
                    coo_type = 'model';
                end
                
                if numel(coo) > N_MAX_POINTS && not(force_label)
                    flag_label = false;
                    flag_label_bg = false;
                end
                
                show_velocities = (not(isempty(id_ref)) || use_displacements) && not(strcmp(arrow_type, 'none'));

                % Manage point coloring
                if ~isempty(fill_val)
                    p_n_col = 201;
                    if isempty(fill_lim)
                        m_m = minMax(fill_val);
                    else
                        m_m = fill_lim;
                    end
                    max_val = max(abs(m_m));
                    if m_m(1) >= 0
                        p_cmap = flipud(Cmap.get('magma', p_n_col));
                        lim = [ m_m(1) max_val];
                    else
                        % map centered to zero
                        p_cmap = Cmap.get('RdBu', p_n_col);
                        lim = [-max_val max_val];
                    end
                    p_color = p_cmap(min(p_n_col,max(1,round((fill_val - lim(1)) / diff(lim) * (p_n_col-1)) + 1)), :);
                end

                % Begin of function ---------------------------------------
                % Get the coordinates of the points to plot
                m_pos = coo.getMedianPos(coo_type);
                [lat, lon, h] = deal(nan(size(m_pos,1),1));
                for p = 1 : numel(m_pos)
                    [lat(p), lon(p), h(p)] = m_pos(p).getGeodetic();
                end
                RED2DEG = 180/pi;
                lat = lat .* RED2DEG;
                lon = lon .* RED2DEG;
                if ~any([lat(:) lon(:)])
                    Logger.getInstance.addWarning('No valid coordinates present, nothing to do');
                else
                    if isempty(fig_handle)
                        if new_fig
                            fh = figure('Visible', 'on');
                        else
                            fh = gcf;
                            hold on;
                        end
                    else
                        fh = fig_handle;
                    end

                    fh_list = fh;
                    if ~isfield(fh.UserData, 'fig_name')
                        fig_name = sprintf('coo_map');
                        fh.UserData.fig_name = fig_name;
                    end

                    Core.getLogger.addMarkedMessage('Preparing coordinates map, please wait...');
                    fh.Color = [1 1 1];

                    [dlat, dlon, dlat_ext, dlon_ext, nwse] = Core_Utils.getMapBoundaries(lat, lon, add_margin);

                    % Manage DTM resolution
                    if strcmp(dtm_resolution, 'auto')
                        patch_dim = sqrt(abs(nwse(1) - nwse(3)) * abs(nwse(2) - nwse(4)));
                        if patch_dim < 1
                            dtm_resolution = '1';
                        elseif patch_dim < 2.5
                            dtm_resolution = '3';
                        else
                            dtm_resolution = 'high';
                        end
                    end
                    setAxis(fh);
                    if flag_no_prj
                         hold on;
                    else
                        ax = setAxis(fh);
                        if new_fig || isempty(ax.Children)
                            switch proj_type
                                case 'mercator'
                                    m_proj('mercator', 'lon', dlon_ext, 'lat', dlat_ext);   % Projection
                                case 'equidistant'
                                    m_proj('equidistant', 'lon', dlon_ext, 'lat', dlat_ext);   % Projection
                                case 'UTM'
                                    m_proj('UTM', 'lon', dlon_ext,'lat', dlat_ext);   % Projection
                                case 'none'
                                otherwise
                                    m_proj('equidistant', 'lon', dlon_ext, 'lat', dlat_ext);   % Projection
                            end
                        end
                    end
                    if show_velocities
                        colorbar;
                    end
                    if new_fig
                        Core_UI.beautifyFig(fh, 'light');
                        drawnow
                    end
                   
                    switch bg_type
                        case {'map'}
                            % White labels in Google maps provides enough contrast
                            if sum(label_color - [0.1 0.1 0.1]) == 0
                                label_color = [1 1 1] -0.001; % Remove 0.001 to trick beautify function
                            end
                            %flag_label_bg = false;

                            if strcmp(proj_type, 'none')
                                if new_fig
                                    addMap('alpha', 0.95, 'lon_lim', dlon_ext, 'lat_lim', dlat_ext, 'provider', map_type);
                                    ylabel('Latitude [deg]');
                                    xlabel('Longitude [deg]');
                                end
                            else
                                img_ggl = [];
                                ax = setAxis(fh);
                                [img_h, lon_ggl, lat_ggl, img_ggl] = addMap('ax', ax, 'alpha', 0.95, 'lon_lim', dlon_ext, 'lat_lim', dlat_ext, 'm_map', true);
                                
                                if contour_lines && max(diff(nwse([3 1])), diff(nwse([2 4]))) < 1 % Add contour
                                    [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', dtm_resolution);
                                    if numel(dtm) > 100
                                        % If the DTM is large enough
                                        dtm = flipud(dtm);
                                        % comment the following line to have bathimetry
                                        dtm(dtm < -50) = nan; % - 1/3 * max(dtm(:));

                                        if any(dtm(:))
                                            step = iif(islogical(contour_lines), iif(diff(minMax(dtm(:))) < 250, 10, iif(diff(minMax(dtm(:))) < 1000, 25, 50)), contour_lines);
                                            if strcmp(proj_type, 'none')
                                                contour(lon_dtm, lat_dtm, dtm)
                                            else
                                                [~, hc] = m_contour(lon_dtm, lat_dtm, dtm, floor(diff(minMax(dtm(:))/step) + 1));
                                            end
                                            hc.LineColor = [1 1 1];
                                            img_h = m_image(lon_ggl, lat_ggl, img_ggl);
                                            img_h.AlphaData = img_h.AlphaData * 0 + 0.6;
                                        end
                                    end
                                end
                            end
                        case {'dtm'}
                            % retrieve external DTM
                            try
                                if sum(label_color - [0.1 0.1 0.1]) == 0
                                    label_color = Core_UI.BLUE; % Remove 0.001 to trick beautify function
                                end
                                flag_label_bg = false;

                                [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', dtm_resolution);
                                dtm = flipud(dtm);
                                % comment the following line to have bathimetry
                                dtm(dtm < 0) = nan; % - 1/3 * max(dtm(:));

                                % uncomment the following line to have colors
                                % colormap(Cmap.adaptiveTerrain(minMax(dtm(:))));

                                setAxis(fh);

                                gray_map = flipud(gray(1000)); colormap(gray_map(150: end, :));
                                drawnow;
                                warning off
                                [shaded_dtm, x, y] = m_shadedrelief(lon_dtm, lat_dtm, dtm, 'nan', [0.98, 0.98, 1], 'gra', 30);
                                warning on
                                %h_dtm = m_pcolor(lon_dtm, lat_dtm, dtm);
                                %h_dtm.CData = shaded_dtm;
                                if strcmp(proj_type, 'none')
                                    img_h = image(lon_dtm, lat_dtm, shaded_dtm); hold on;
                                    if flag_new
                                        fh.Children(end).YDir = 'normal';
                                        axis equal; axis tight;
                                        xlim(dlon);
                                        ylim(dlat);
                                    end
                                else
                                    img_h = m_image(lon_dtm, lat_dtm, shaded_dtm); hold on;
                                end
                                if contour_lines && any(dtm(:)) && numel(dtm) > 100 % Add contour
                                    img_h.AlphaData = img_h.AlphaData .* 0.70; % 70% of alpha on DTM
                                    step = iif(islogical(contour_lines), iif(diff(minMax(dtm(:))) < 250, 10, 50), contour_lines);
                                    if strcmp(proj_type, 'none')
                                        contour(lon_dtm, lat_dtm, dtm)
                                    else
                                        [~, hc] = m_contour(lon_dtm, lat_dtm, dtm, floor(diff(minMax(dtm(:))/step) + 1));
                                    end
                                    hc.LineColor = 0.2 + [0 0 0];
                                end
                            catch ex
                                % use ETOPO1 instead
                                try
                                    m_etopo2('shadedrelief','gradient', 3);
                                catch
                                    bg_type = none;
                                end
                            end
                        otherwise
                            if contour_lines % Add contour
                                [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', dtm_resolution);
                                if numel(dtm) > 100
                                    % If the DTM is large enough
                                    dtm = flipud(dtm);
                                    % comment the following line to have bathimetry
                                    dtm(dtm < -50) = nan; % - 1/3 * max(dtm(:));

                                    if any(dtm(:))
                                        step = iif(islogical(contour_lines), iif(diff(minMax(dtm(:))) < 250, 10, iif(diff(minMax(dtm(:))) < 1000, 25, 50)), contour_lines);
                                        [~, hc] = m_contour(lon_dtm, lat_dtm, dtm, floor(diff(minMax(dtm(:))/step) + 1));
                                        hc.LineColor = [0.5 0.5 0.5];
                                    end
                                end
                            end
                    end

                    % read shapefile if requested
                    if (~strcmp(shape,'none'))
                        if (~strcmp(shape,'coast')) && (~strcmp(shape,'fill'))
                            if (strcmp(shape,'10m'))
                                M = m_shaperead('countries_10m');
                            elseif (strcmp(shape,'30m'))
                                M = m_shaperead('countries_30m');
                            else
                                M = m_shaperead('countries_50m');
                            end
                            [x_min, y_min] = m_ll2xy(min(lon_lim), min(lat_lim));
                            [x_max, y_max] = m_ll2xy(max(lon_lim), max(lat_lim));
                            for k = 1 : length(M.ncst)
                                lam_c = M.ncst{k}(:,1);
                                ids = lam_c <  min(lon);
                                lam_c(ids) = lam_c(ids) + 360;
                                phi_c = M.ncst{k}(:,2);
                                [x, y] = m_ll2xy(lam_c, phi_c);
                                if sum(~isnan(x))>1
                                    x(find(abs(diff(x)) >= abs(x_max - x_min) * 0.90) + 1) = nan; % Remove lines that occupy more than th 90% of the plot
                                    line(x,y,'color', [0.3 0.3 0.3]);
                                end
                            end
                        else
                            if (strcmp(shape,'coast'))
                                m_coast('line','color', lineCol);
                            else
                                m_coast('patch',lineCol);
                            end
                        end
                    end

                    hold on;
                    if  ~strcmp(proj_type, 'none')
                        m_grid('box','fancy','tickdir','in', 'fontsize', Core_UI.getFontSize(11));
                        % m_ruler(1.1, [.05 .40], 'tickdir','out','ticklen',[.007 .007], 'fontsize',14);
                        m_ruler([.7 1], -0.075, 'tickdir','out','ticklen',[.007 .007], 'fontsize',Core_UI.getFontSize(9));
                    else
                        fh.Children(end).XAxis.TickLength = [0.005, 0.005];
                        fh.Children(end).YAxis.TickLength = [0.005, 0.005];
                        grid(fh.Children(end), 'on')
                        %grid(fh.Children(end), 'minor')
                    end
                    drawnow;

                    % Display points
                    if strcmp(proj_type, 'none')
                        x = lon;
                        y = lat;
                    else
                        [x, y] = m_ll2xy(lon, lat);
                    end

                    if isfield(fh.UserData, 'x')
                        fh.UserData.x = [fh.UserData.x(:); x(:)];
                        fh.UserData.y = [fh.UserData.y(:); y(:)];
                    else
                        fh.UserData.x = x(:);
                        fh.UserData.y = y(:);
                    end
                    % Draw labels bg (behind the points ---------------------------

                    if flag_label && flag_label_bg
                        % Label BG (in background w.r.t. the point)
                        for r = 1 : numel(coo)
                            name = upper(coo(r).name);
                            name = pad(name(1:min(4,numel(name))), 4, 'both', ' ');
                            hlbl_bg(r) = text(x(r), y(r), ['  ' name '  '], ...
                                'FontWeight', 'bold', 'FontSize', label_size, 'Color', [1 1 1], ...
                                'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                                'Margin', 2, 'LineWidth', 1, ...
                                'HorizontalAlignment','left');
                        end
                    end

                    % Plot velocities
                    if show_velocities
                        % Recompute velocities with the current reference ---------
                        v_enu = zeros(numel(coo), 3);
                        if animate
                            v_enu_final = zeros(numel(coo), 3);
                        end
                        coo = coo.getCopy;
                        if not(isempty(id_ref))
                            coo.setNewRef(id_ref);
                        end
                        if use_displacements
                            if isempty(id_ref)
                                idr = [];
                            else
                                idr = id_ref;
                            end
                            % any coordinate is good as local ref
                            for r = 1 : numel(coo)
                                if isempty(idr)
                                    enu = coo(r).getLocal(coo(r).getMedianPos, coo_type);
                                    %enu = coo(r).getLocal(coo(r).getElement(1), coo_type);
                                else
                                    enu = coo(r).getLocal(coo(idr).getMedianPos, coo_type);
                                    %enu = coo(r).getLocal(coo(idr).getElement(1), coo_type);
                                end
                                enu = enu(not(isnan(enu(:,1))), :); % keep only valid values
                                enu = bsxfun(@minus, enu, enu(1,:));% displacement relative to the first valid epoch
                                if animate
                                    % get the maximum displacement %the sign will not be important
                                    v_enu_final(r,:) = enu(end,:);
                                    v_enu(r,:) = [max(abs(minMax(enu(:,1)))), max(abs(minMax(enu(:,2)))), max(abs(minMax(enu(:,3))))];
                                else
                                    v_enu(r,:) = enu(end,:);
                                end
                            end
                        else
                            for r = 1 : numel(coo)
                                if (r ~= id_ref)
                                    if recompute_vel
                                        coo(r).v_xyz = [];
                                        [~, v_enu(r,:)] = coo(r).getLocal(coo(id_ref).getMedianPos, coo_type);
                                    else
                                        tmp_vel = mean(coo(r).v_xyz, 1);
                                        if isempty(tmp_vel)
                                            tmp_vel = [0 0 0];
                                        end
                                        v_enu(r,:) = tmp_vel;
                                        clear tmp_vel
                                    end
                                end
                            end
                        end
                        flag_speed = not(use_displacements); % the user wants to visualize displacements instead of speed
                        if any(v_enu(:))
                            [ah, max_size, lim_vel] = Coordinates.insertArrowsLegend(fh, v_enu, vel_lim, arrow_type, flag_speed);
                            if animate
                                ah_list = cell(numel(x),1);
                                v_enu = v_enu_final;
                                ph_list = [];
                            else
                                [ah_list, ph_list] = Coordinates.insertArrows(fh, id_ref, x, y, v_enu, max_size, lim_vel, cvel_lim, cmap, arrow_type, true);
                            end
                            colormap(cmap); caxis([0 cvel_lim(1 + (strcmp(arrow_type, 'up')))]); cb = colorbar; title(cb, '[mm]', 'FontSize', 8);
                        else
                            % No velocities, no arrows!
                            arrow_type = 'none';
                            animate = false;
                        end
                    
                    % Set points style --------------------------------------------

                        for r = 1 : numel(coo)
                            if not(isempty(ph_list)) && not(isempty(ph_list{r}))
                                ph(r) = ph_list{r};
                                ph(r).MarkerSize = point_size;
                                ph(r).LineWidth = 2;
                                if ismember(r, id_ref)
                                    ph(r).Marker = '^';
                                    ph(r).MarkerEdgeColor = max(0, Core_UI.ORANGE - 0.1);
                                    ph(r).MarkerFaceColor = min(1, Core_UI.ORANGE + 0.1);
                                else
                                    ph(r).Marker = 'o';
                                    ph(r).MarkerEdgeColor = max(0, point_color - 0.1);
                                    ph(r).MarkerFaceColor = min(1, point_color + 0.1);
                                end
                            end
                        end
                    else
                        % Draw points -------------------------------------------------

                        if flag_tooltip
                            for r = 1 : numel(coo)
                                if ismember(r, id_ref)
                                    if ~isempty(fill_val)
                                        ph(r) = plot(x(r), y(r), '^', 'MarkerSize', point_size, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', p_color(r,:), 'LineWidth', 1); hold on;
                                    else
                                        ph(r) = plot(x(r), y(r), '^', 'MarkerSize', point_size, 'MarkerEdgeColor', max(0, Core_UI.ORANGE - 0.1), 'MarkerFaceColor', min(1, Core_UI.ORANGE + 0.1), 'LineWidth', 2); hold on;
                                    end
                                else
                                    if ~isempty(fill_val)
                                        ph(r) = plot(x(r), y(r), 'o', 'MarkerSize', point_size, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', p_color(r,:), 'LineWidth', 1); hold on;
                                    else
                                        ph(r) = plot(x(r), y(r), 'o', 'MarkerSize', point_size, 'MarkerEdgeColor', max(0, point_color - 0.1), 'MarkerFaceColor', min(1, point_color + 0.1), 'LineWidth', 2); hold on;
                                    end
                                end
                            end
                        else
                            p_list = 1 : numel(coo);
                            r = setdiff(p_list, id_ref);
                            if any(r)
                                if ~isempty(fill_val)
                                    ph(1) = plot(x(r), y(r), 'o', 'MarkerSize', point_size, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', p_color(r,:), 'LineWidth', 1); hold on;
                                else
                                    ph(1) = plot(x(r), y(r), 'o', 'MarkerSize', point_size, 'MarkerEdgeColor', max(0, point_color - 0.1), 'MarkerFaceColor', min(1, point_color + 0.1), 'LineWidth', 2); hold on;
                                end
                            end
                            r = intersect(p_list, id_ref);
                            if any(r)
                                if ~isempty(fill_val)
                                    ph(2) = plot(x(r), y(r), '^', 'MarkerSize', point_size, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', p_color(r,:), 'LineWidth', 1); hold on;
                                else
                                    ph(2) = plot(x(r), y(r), '^', 'MarkerSize', point_size, 'MarkerEdgeColor', max(0, Core_UI.ORANGE - 0.1), 'MarkerFaceColor', min(1, Core_UI.ORANGE + 0.1), 'LineWidth', 2); hold on;
                                end
                            end
                        end
                    end                    

                    % Draw labels -------------------------------------------------

                    if flag_label
                        lbl_pos = nan(numel(coo), 4);
                        for r = 1 : numel(coo)
                            name = upper(coo(r).name);
                            name = pad(name(1:min(4,numel(name))), 4, 'both', ' ');
                            hlbl(r) = text(x(r), y(r), ['  ' name '  '], ...
                                'FontWeight', 'bold', 'FontSize', label_size, 'Color', label_color, ...
                                'Margin', 2, 'LineWidth', 2, ...
                                'HorizontalAlignment','left');
                            lbl_pos(r,:) = [hlbl(r).Position(1:2) hlbl(r).Extent(3:4)];
                            if not(flag_label_bg)
                                hlbl(r).UserData = struct('keep_color', true);
                            end
                        end

                        % Check lbl overposition
                        if numel(lbl_pos(:,1)) < 90
                            if flag_fix_label
                                if flag_label_bg
                                    Coordinates.fixLabels(lbl_pos, [hlbl(:) hlbl_bg(:)]);
                                else
                                    Coordinates.fixLabels(lbl_pos, hlbl(:));
                                end
                            end
                        end
                        setAxis(fh);
                    end

                    % Draw infoboxes ----------------------------------------------

                    if flag_tooltip
                        ax = setAxis(fh);
                        old_state = ax.Units;
                        ax.Units = 'points';
                        scale_x = diff(xlim) / ax.Position(3);
                        if use_displacements
                            unit = 'mm';
                            mov_title = 'Displacements';
                        else
                            unit = 'mm/y';
                            mov_title = 'Velocities';
                        end
                        for r = 1 : numel(coo)
                            if ismember(r, id_ref)
                                [lat, lon, h_ellips, h_ortho] = coo(r).getMedianPos.getGeodetic();
                                str = {sprintf('Master station: %s', coo(r).name), ...
                                    iif(length(coo(r).description) <= 4, sprintf('%s', coo(r).name_v3), sprintf('%s', coo(r).description)), ...
                                    sprintf('Lat: %.5f deg', lat / pi * 180), ...
                                    sprintf('Lon: %.5f deg', lon / pi * 180), ...
                                    sprintf('El: %.2f m', h_ellips), ...
                                    sprintf('El: %.2f m (orthometric)', h_ortho)};
                            else
                                [lat, lon, h_ellips, h_ortho] = coo(r).getMedianPos.getGeodetic();
                                str = {sprintf('Station: %s', coo(r).name), ...
                                    iif(length(coo(r).description) <= 4, sprintf('%s', coo(r).name_v3), sprintf('%s', coo(r).description)), ...
                                    sprintf('Lat: %.5f deg', lat / pi * 180), ...
                                    sprintf('Lon: %.5f deg', lon / pi * 180), ...
                                    sprintf('El: %.2f m', h_ellips), ...
                                    sprintf('El: %.2f m (orthometric)', h_ortho)};
                                if not(isempty(id_ref))
                                    loc_enu = coo(r).getMedianPos.getLocal(coo(id_ref).getMedianPos);
                                    bsl = hypot(loc_enu(1), loc_enu(2));
                                    if bsl > 1000
                                        bsl = sprintf('%.3f km', bsl * 1e-3);
                                    else
                                        bsl = sprintf('%.3f m', bsl);
                                    end
                                    str = [str, {'--------------------------------', ...
                                        sprintf('Baseline: %s - %s', coo(r).name, coo(id_ref).name), ...
                                        sprintf('Len: %s (median)', bsl), ...
                                        sprintf('El. diff: %.3f m (median)', loc_enu(3))}];
                                end

                                if show_velocities
                                    % if velocities have been computed
                                    str = [str {'', sprintf('%s %s', coo(r).name, mov_title), ...
                                        sprintf('Planar: %.1f %s', sqrt(v_enu(r,1).^2 + v_enu(r,2).^2) .* 1e3, unit), ...
                                        sprintf('Nord: %.1f %s', v_enu(r,2) .* 1e3, unit), ...
                                        sprintf('East: %.1f %s', v_enu(r,1) .* 1e3, unit), ...
                                        sprintf('Up: %.1f %s', v_enu(r,3) .* 1e3, unit)}];
                                end
                            end
                            if ~isempty(fill_val)
                                str = [str {sprintf('Value: %g', fill_val(r))}];
                            end
                            if verLessThan('matlab', '99.9') % (this is faster)
                                txt = text(x(r), y(r), str, ...
                                    'BackgroundColor', [1 1 1], 'EdgeColor', [0.3 0.3 0.3], ...
                                    'Margin', 5, ...
                                    'FontSize', 9, ...
                                    'Visible', 'off', ...
                                    'Interpreter', 'none');
                                ph(r).UserData = struct('label', txt);
                            else
                                % Add datatip (this is too slow)
                                n_el = numel(str);
                                dt = datatip(ph(r));
                                ph(r).DataTipTemplate.DataTipRows = dataTipTextRow('NONE','');
                                for i = 1:n_el
                                    ph(r).DataTipTemplate.DataTipRows(i).Label = str{i};
                                end
                                dt.Visible = 'off';
                                ph(r).UserData = struct('datatip', dt);
                            end
                        end
                        if isfield(fh.UserData, 'ph')
                            fh.UserData.ph = [fh.UserData.ph(:); ph(:)];
                        else
                            fh.UserData.ph = ph(:);
                        end
                        ax = setAxis(fh);
                        fh.WindowButtonMotionFcn = {@Coordinates.mouseMove, ax, [fh.UserData.x(:) fh.UserData.y(:)], fh.UserData.ph(:)};
                        ax.Units = old_state;
                    end

                    % Complete fig presentation -----------------------------------
                    str_title = 'Map of GNSS stations';
                    th = title(sprintf('%s\\fontsize{5} \n', str_title));
                    if new_fig
                        Core_UI.addExportMenu(fh); Core_UI.addBeautifyMenu(fh); Core_UI.beautifyFig(fh, 'light');
                        fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                    end
                                        
                    if new_fig
                        if ~isempty(fill_val)
                            colormap(p_cmap);
                            caxis(lim);
                            colorbar;
                        end
                    end

                    % Animation ---------------------------------------------------

                    if animate && use_displacements
                        times = [coo.time];
                        [sync_time, id_sync] = times.getSyncedTime();
                        enu = nan(sync_time.length, 3, numel(coo));
                        for r = 1 : numel(coo)
                            % any coordinate is good as local ref
                            tmp = coo(r).getLocal(coo(idr).getMedianPos, coo_type);
                            % rem duplicates
                            tmp(find(diff(coo(r).time.getMatlabTime) == 0), :) = [];
                            enu(not(isnan(id_sync(:,r))), :, r) = tmp;
                            id_first = find(not(isnan(enu(:, 1, r))), 1, 'first');

                            enu(:, :, r) = bsxfun(@minus, enu(:, :, r), enu(id_first, :, r));
                            %enu(:, :, r) = -bsxfun(@minus, enu(:, :, r), mean(enu(:, :, r), 'omitnan'));
                        end

                        % Animate arrows
                        max_fps = 10; % max number of fps
                        max_t = 30;   % 10 seconds limit
                        max_n_frames = sync_time.length;

                        % First export
                        ani_export = ani_export;
                        if ani_export(1) == 'g' || ani_export(1) == 'v' % if export mode is gif or video
                            th.String = sprintf('%s\n%s\\fontsize{5} \n', str_title, sync_time.getEpoch(1).toString('yyyy-mm-dd HH:MM'));
                            drawnow
                            if not(isempty(ani_file_name))
                                file_name = ani_file_name;
                            else
                                file_name = fullfile(Core.getState.getOutDir, ['Ani' fig_name]);
                                if ani_export(1) == 'g'
                                    file_name = [file_name '.gif'];
                                end
                            end

                            ss_image = 2; %% sub_sample image
                            frameset = {};
                            frame = getframe(fh);
                            frame = frame(1:ss_image:end, 1:ss_image:end,:);
                            if ani_export(1) == 'g'
                                im = frame2im(frame);
                                [imind, cm] = rgb2ind(im,256);
                                % Write to the GIF File
                                imwrite(imind, cm, file_name, 'gif', 'Loopcount', inf, 'DelayTime', 1/max_fps);
                            else
                                frameset{1} = frame; % subsample (1:2)
                            end
                            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                            Core.getLogger.addMarkedMessage('Exporting video');
                            fprintf('%5d/%5d', 0, 99999);
                        end

                        % Compute rate to stay under 10s and 25 fps
                        ss_rate = ceil(max_n_frames / max_fps / max_t);
                        epoch_list = (ss_rate + 1) : ss_rate : max_n_frames;
                        v_enu = zeros(numel(coo), 3);
                        for i = 1 : numel(epoch_list)
                            e = epoch_list(i);
                            th.String = sprintf('%s\n%s\\fontsize{5} \n', str_title, sync_time.getEpoch(e).toString('yyyy-mm-dd HH:MM'));
                            for r = 1 : numel(coo)
                                for h = ah_list{r}
                                    delete(h);
                                end
                                if any(enu(e, :, r))
                                    v_enu(r,:) = enu(e, :, r);
                                end
                            end

                            ah_list = Coordinates.insertArrows(fh, id_ref, x, y, v_enu, max_size, lim_vel, cvel_lim, cmap, arrow_type, flag_speed);
                            drawnow

                            if ani_export(1) == 'g' || ani_export(1) == 'v'
                                fprintf('%s%5d/%5d',char(8 * ones(1,11)), i, numel(epoch_list));
                                frame = getframe(fh);
                                frame = frame(1:ss_image:end, 1:ss_image:end,:);
                                if ani_export(1) == 'g'
                                    im = frame2im(frame);
                                    [imind, cm] = rgb2ind(im, 256);
                                    % Write to the GIF File
                                    if i == numel(epoch_list)
                                        imwrite(imind, cm, file_name, 'gif', 'WriteMode', 'append', 'DelayTime', 5);
                                    else
                                        imwrite(imind, cm, file_name, 'gif', 'WriteMode', 'append', 'DelayTime', 1/max_fps);
                                    end
                                else
                                    frameset{i + 1} = frame; % subsample (1:2)
                                end
                            end
                        end

                        % Finish file save
                        if ani_export(1) == 'g'
                            Core.getLogger.addStatusOk(sprintf('"%s" done ^_^', fullfile(Core.getState.getOutDir, file_name)));
                        elseif ani_export(1) == 'v'
                            fprintf('%s',char(8 * ones(1,11)));
                            if ismac() || ispc()
                                % Better compression on Mac > 10.7 and Win > 7
                                video_out = VideoWriter(fullfile(Core.getState.getOutDir, ['AniMap' fig_name '.mp4']), 'MPEG-4');
                            else
                                % Linux doesn't have mp4 compression avaiable
                                video_out = VideoWriter(fullfile(Core.getState.getOutDir, ['AniMap' fig_name '.avi']));
                            end
                            video_out = VideoWriter(file_name);

                            video_out.FrameRate = max_fps;
                            video_out.Quality = 91;
                            open(video_out);
                            for i = 1 : numel(frameset)
                                writeVideo(video_out, frameset{i});
                            end
                            close(video_out);
                            Core.getLogger.addStatusOk(sprintf('"%s" done ^_^', fullfile(Core.getState.getOutDir, video_out.Filename)));
                        end

                    end
                end
            end
        end
        
        function [smooth_rh, time_out, err] = getReflectorHeightInterp(coo, spline_base, rate, max_merge)
            % Compute an interpolation of the rreflector height (useful in the cases of raw sample saved)
            %
            % INPUT
            %   coo           a single coordinate
            %   spline_base   default 3 hours (in seconds)
            %   out_rate      default 15 minutes (in seconds)
            %   max_merge     apply output mask (max n epochs for extrapolation)
            %
            % OUTPUT
            %   smooth_rh     spline interpolated reflector height
            %   time_out      time of the interpolated points
            %   err           approximate estimated error as interpolated dispersion of the data
            %
            % SYNTAX
            %   [smooth_rh, time_out, err] = getReflectorHeightInterp(coo, spline_base, rate)

            if nargin < 2 || isempty(spline_base)
                spline_base = 86400/24*3;
            end
            if nargin < 3 || isempty(rate)
                rate = 900;
            end
            if nargin < 4 || isempty(max_merge)
                max_merge = 8; % 2 hours
            end

            % Get the data
            rh = coo.reflectometry;
            % Compute the output time
            ref_time = floor(rh.time.first.getMatlabTime);
            time_in = rh.time.getRefTime(ref_time);
            start_time = floor(time_in(1)/rate)*rate;
            stop_time =  ceil(time_in(end)/rate)*rate;
            time_out = (start_time:rate:stop_time)';
            % prepare data
            data_filt = double([rh.value max(0.05, rh.std).^2]);
            
            % compute a mask
            mask_out = false(size(time_out));
            mask_out(round(time_in/rate)+1 - floor(time_in(1)/rate)) = true;
            mask_out = flagMerge(mask_out, max_merge);

            % Compute a sort of std, based on the dispersion of the observations
            t_span = time_out(end)-time_in(1);
            n = ceil(t_span / spline_base - 1);
            t_dispersion = linspace(time_out(1), time_out(end), n);
            step = t_span/(n-1);
            idd = min(round(time_in/step)+1, numel(time_out));
            err = nan(n,1);
            for i = 1:n
                err(i) = std(rh.value(idd == i));
            end
            id_ok = ~isnan(err);
            err = interp1(t_dispersion(id_ok), err(id_ok), time_out, 'pchip'); % piece wise cubic interpolation

            % compute interpolation
            m = robAdj(data_filt(:,1)');
            
            % compute splines for regularization
            s = 10*robStd(data_filt(:,1))/2;
            % calibrate them on the biggest gap, min 4 spline_bas, max 1 day
            max_gap = min(86400, max(4*spline_base, max(diff(find(mask_out))) * rate));
            [~, ~, ~, smooth_rh] = splinerMat(time_in, [data_filt(:,1) - m data_filt(:,2)], max_gap, 1, time_out);
            % keep the regularization data as values in the gaps
            % use s of error for them
            data_filt = [data_filt; [smooth_rh(~mask_out) + m ones(sum(~mask_out),1) * s^2]];
            time_in = [time_in; time_out(~mask_out)];
            % sort the input with added pseudo obs
            [time_in, ids] = sort(time_in);
            data_filt = data_filt(ids,:);
            [~, ~, ~, smooth_rh] = splinerMat(time_in, [data_filt(:,1) - m data_filt(:,2)], spline_base, 0.5, time_out);
            smooth_rh = smooth_rh + m;
            % mask values in the gaps
            smooth_rh(~mask_out) = nan;
            err(~mask_out) = nan;

            % Convert back time_out to GPS_Time
            time_out = GPS_Time(ref_time + time_out(:)/86400);
            
            % check for divergence
            % if the interpolation is less than the min value or higher than the higher value by more the 1/2 of
            % the observation variation, mark it as outlier
            data_lim = minMax(data_filt(:,1));
            tollerance = diff(data_lim)/2;
            smooth_rh(smooth_rh < data_lim(1) - tollerance | smooth_rh > data_lim(2) + tollerance) = nan;
        end

        function fh = showReflectorHeightFromRaw(coo, flag_height)
            % Show the results from raw reflectometry data as height
            %
            % INPUT
            %   coo         Coordinate object containing raw reflectometry data
            %   flag_height if true, display height; otherwise, display distance from the antenna
            %
            % SYNTAX
            %   fh = coo.showReflectorHeightFromRaw(flag_height)

            % Check input arguments
            if nargin < 2 || isempty(flag_height)
                flag_height = true;
            end

            % Get interpolated reflector height data on a regular grid
            [rh_value, time_data, err] = coo.getReflectorHeightInterp(86400/12, 900);
            coo_name = coo.getNameV3;

            % Get original reflectometry time and value
            raw_time = coo.reflectometry.time;
            raw_value = coo.reflectometry.value;
            [~, ~, h_ell, h_ort] = coo.getGeodetic;
            time_coo = coo.getTime;
                    
            % Adjust original values if flag_height is true
            if flag_height
                if time_coo.length > 1
                    raw_value = interp1(time_coo.getMatlabTime, h_ort, raw_time.getMatlabTime, 'spline') - raw_value;
                    rh_value = interp1(time_coo.getMatlabTime, h_ort, time_data.getMatlabTime, 'spline') - rh_value;
                else
                    raw_value = h_ort - raw_value;
                    rh_value = h_ort - rh_value;
                end
            end

            % Create a new figure
            fh = figure('Visible', 'off');
            fig_name = sprintf('%s RH', coo_name);
            fh = figure('Visible', 'off');  Core_UI.beautifyFig(fh);
            fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
            if flag_height
                fig_name = sprintf('IR_Results_%s_%s-%s', coo_name, time_data.first.toString('yyyymmdd_HHMM'), time_data.last.toString('yyyymmdd_HHMM'));
            else
                fig_name = sprintf('IR_Results_distance_%s_%s-%s', coo_name, time_data.first.toString('yyyymmdd_HHMM'), time_data.last.toString('yyyymmdd_HHMM'));
            end
            fh.UserData = struct('fig_name', fig_name);
            set(fh,'defaultAxesColorOrder',Cmap.getColor([15 70],100, 'viridis'));

            Core_UI.beautifyFig(fh);

            % Get the datetime array for x-axis from the interpolated data
            t = time_data.getDateTime();

            % Add the patch for standard deviation as a confidence interval
            id_ok = true(size(rh_value));
            Core_Utils.plotConfBand(t, rh_value, 3*err, Core_UI.BLUE, 0.1);
            Core_Utils.plotConfBand(t, rh_value, err, Core_UI.BLUE, 0.1);
            fh.Visible = 'on';
            hold on;

            % Plot original observations as light gray dots
            t_raw = raw_time.getDateTime();
            plot(t_raw, raw_value, '.', 'Color', ones(1,3)*0.7, 'DisplayName', 'RAW Observations');
            
            % Plot interpolated data
            plot(time_data.getDateTime(), rh_value, '.-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Interpolated Water Level', 'Color', Core_UI.BLUE);

            % Set labels and title based on flag_height
            if flag_height
                ylabel('Reflector Height [m]');
                title(sprintf('GNSS-IR results: %s\\fontsize{5}\n', coo_name));
            else
                ylabel('Reflector Distance [m]');
                title(sprintf('GNSS-IR results: %s\\fontsize{5}\n', coo_name));
                set(gca,'Ydir','reverse'); % Reverse Y-axis if displaying distance
            end

            % Add a legend
            legend('Location', 'best');

            % Additional plotting setup
            grid on;
            xlim(minMax([t; t_raw])); % Set the x-axis limits to include both original and interpolated data

            % Beautify and export options for the figure
            Core_UI.addBeautifyMenu(fh);
            Core_UI.addExportMenu(fh);
            Core_UI.beautifyFig(fh);
            fh.Visible = iif(Core_UI.isHideFig, 'off', 'on');

            drawnow;
        end


        function fh = showReflectorHeight(coo, flag_height, flag_info)
            % Show the results from reflectometry height (default) | distance
            %
            % INPUT
            %   coo         if it is a list shows the mean
            %   flag_height if true display heigtht otherwise distance from the antenna
            %   flag_info   works in single coo mode, add n_data subplot
            %
            % SYNTAX
            %   fh = coo.showReflectorHeight(flag_height, flag_info)
            fh = [];
            min_n_obs = 8;
            if nargin < 2 || isempty(flag_height)
                flag_height = true;
            end
            if numel(coo) > 1
                % With multiple coo compute mean
                [rh_data, time_data, data_sync, n_obs] = coo.getMeanReflectionHeight(flag_height);
                if all(n_obs == 1)
                    min_n_obs = 1;
                end
                rh_data(all(n_obs < min_n_obs,2),:) = nan;
                coo_name = 'mean reflector height [masl]';
                flag_ok = true;
                flag_info = false;
            else
                % Single coo mode
                flag_ok = not(isempty(coo.reflectometry.value));
                if flag_ok
                    [~, ~, h_ell, h_ort] = coo.getGeodetic;
                    time_coo = coo.getTime;
                    time_data = coo.getReflectometryTime;
                    if flag_height
                        if time_coo.length > 1
                            rh_data = interp1(time_coo.getMatlabTime, h_ort, time_data.getMatlabTime, 'spline') - coo.reflectometry.value;
                        else
                            rh_data = h_ort - coo.reflectometry.value;
                        end
                    else
                        rh_data = coo.reflectometry.value;
                    end
                    if any(coo.reflectometry.n_obs(:) > 1)
                        rh_data(coo.reflectometry.n_obs(:) < min_n_obs) = nan;
                    end
                    coo_name = coo.getNameV3;
                end
                if nargin < 3
                    flag_info = true;
                end
            end
            if flag_ok
                fig_name = sprintf('%s RH', coo_name);
                fh = figure('Visible', 'off');  Core_UI.beautifyFig(fh);
                fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                if flag_height
                    fig_name = sprintf('IR_Results_%s_%s-%s', coo_name, time_data.first.toString('yyyymmdd_HHMM'), time_data.last.toString('yyyymmdd_HHMM'));
                else
                    fig_name = sprintf('IR_Results_distance_%s_%s-%s', coo_name, time_data.first.toString('yyyymmdd_HHMM'), time_data.last.toString('yyyymmdd_HHMM'));
                end
                fh.UserData = struct('fig_name', fig_name);
                set(fh,'defaultAxesColorOrder',Cmap.getColor([15 70],100, 'viridis'));
                if flag_info
                    ax = subplot(4,1,1:3);
                end
                t = time_data.getDateTime();
                if numel(coo) > 1
                    data_sync(n_obs < 5) = nan; % If the reflector height have been computed with less than 5 values ignore it
                    for i = 1:size(data_sync,2)
                        plotSep(t, data_sync(:,i), '.-', 'LineWidth', 1, 'MarkerSize', 2, 'Color', [0.5 0.5 0.5]); hold on;
                    end
                end
                plotSep(t, rh_data, '.-', 'LineWidth', 2, 'MarkerSize', 8);
                if not(flag_info)
                    ax = setAxis(fh);
                end
                %setTimeTicks();
                grid on;
                if flag_height
                    ylabel('Reflector Height [m]');
                else
                    ylabel('Reflector Distance [m]');
                end
                title(sprintf('GNSS-IR results: %s\\fontsize{5}\n', coo_name));
                if numel(coo) == 1
                    ax2 = setAxis(fh);
                    yyaxis right;
                    ylabel('Obs results std [cm]');
                    Core_Utils.patchSep(t, coo.reflectometry.std(:) .* 1e2, Core_UI.getColor(4), 'FaceColor', Cmap.getColor(70,100, 'viridis'), 'EdgeColor','none','FaceAlpha',0.3,'HandleVisibility','off'); hold on;
                    yl = max([0 perc(coo.reflectometry.std(:) .* 1e2, 0.997)], 2);
                    if logical(yl(end))
                        ylim(ax2, [0 yl(2)*3/2]);
                    end
                    yyaxis left
                end

                if ~flag_height
                    set(ax,'Ydir','reverse')
                end
                xlim(minMax(t));
                
                if flag_info
                    ax3 = subplot(4,1,4);
                    plot(t, coo.reflectometry.n_obs, '.--', 'LineWidth', 2, 'MarkerSize', 5, 'Color', Cmap.getColor(10,100, 'viridis'));
                    ylabel('# obs');
                    yl = ylim();
                    grid on;
                    ylim([0 yl(2)]);
                    xlim(minMax(t));
                    linkaxes([ax ax2 ax3], 'x');
                end
                Core_UI.addBeautifyMenu(fh);
                Core_UI.addExportMenu(fh);
                Core_UI.beautifyFig(fh);
                fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
            end
        end
    end
    
    % =========================================================================
    %%    STATIC AUXILLIARY
    % =========================================================================
    
    methods (Static, Access = 'private')
        function status = safeSave(out_file_name, coo)
            if ~(exist(out_file_name, 'file') == 2)
                % If the coordinate file not exist
                [base, fname, fext] = fileparts(out_file_name);
                if strcmp(fext,'.mat')
                    out_file_name = fullfile(base, fname);
                end
                % check if a backup exist
                if (exist([out_file_name '.bk.mat'], 'file') == 2)
                    % promote the backup to primary file
                    [status, msg] = movefile([out_file_name '.bk.mat'], [out_file_name '.mat']);
                end
            end
            % try to see if the original file is not corrupted
            if (exist(out_file_name, 'file') == 2)
                try
                    tmp = load(out_file_name, 'coo'); clear tmp;
                    [base, fname, fext] = fileparts(out_file_name);
                    if strcmp(fext,'.mat')
                        out_file_name = fullfile(base, fname);
                    end
                    [status, msg] = movefile([out_file_name '.mat'], [out_file_name '.bk.mat'], 'f');
                    if ~status
                        Core.getLogger.addError(sprintf('Renaming "%s" failed: %s', out_file_name, msg));
                    end
                catch ex
                    Core.getLogger.addError(sprintf('Renaming "%s" failed: %s', out_file_name, ex.message));
                end
            end
            try
                save(out_file_name, 'coo'); clear tmp;
                tmp = load(out_file_name, 'coo');
                if (exist([out_file_name '.bk.mat'], 'file') == 2)
                    delete([out_file_name '.bk.mat']);
                end
            catch ex
                Core.getLogger.addError(sprintf('Saving "%s" failed: %s', out_file_name, ex.message));
                % check if a backup exist and restore it
                if (exist([out_file_name '.bk.mat'], 'file') == 2)
                    % promote the backup to primary file
                    [status, msg] = movefile([out_file_name '.bk.mat'], out_file_name);
                end
            end
        end

        function [ah, max_size, lim_vel] = insertArrowsLegend(fh, v_enu, vel_lim, arrow_type, flag_speed)
            % Insert velocity arrows in a map  of stations
            %
            % SYNTAX
            %   [ah, max_size, lim_vel] = Coordinates.insertArrowsLegend(fh, v_enu, vel_lim, arrow_type, flag_speed)
            
            if nargin < 5 || isempty(flag_speed)
                flag_speed = true;
            end
            h_scale = 15;
            ax = setAxis(fh);
            max_v = max(abs(v_enu));
            
            % if there are moving points
            xl = xlim();
            yl = ylim();
            old_state = ax.Units;
            ax.Units = 'points';
            scale_x = diff(xlim) / ax.Position(3); % points to X-axis coordinates
            scale_y = diff(ylim) / ax.Position(4); % points to Y-axis coordinates
            
            max_size_x =  ax.Position(3) ./ 5 * scale_x; % 2 / 5 of the map in map size
            max_size_y =  ax.Position(4) ./ 5 * scale_y; % 2 / 5 of the map in map size
            lim_vel_x = 0;
            
            if any(max_v)
                if (strcmp(arrow_type, 'up'))
                    lim_vel = max(vel_lim(2,2), (min(vel_lim(1,2), max_v(3) * 1e3)));
                else
                    lim_vel_x = min(vel_lim(1,1), (max(vel_lim(1,2), max_v(1) * 1e3)));
                    lim_vel_y = min(vel_lim(2,1), (max(vel_lim(1,2), max_v(2) * 1e3)));
                    lim_vel = max(lim_vel_x, lim_vel_y);
                end
                % Draw the arrow "Legend" -----------------------------
                
                v_legend = max(5,floor((lim_vel / max_size_x * 80 * scale_x)/20)*20);
                
                hr = rectangle('Position', [xl(1) + 20 * scale_x, yl(1) + 20 * scale_y, 20 * scale_y + v_legend / lim_vel * max_size_x, 30 * scale_y], 'FaceColor', Core_UI.WHITE);
                if flag_speed
                    ht = text(xl(1) + 30 * scale_x, yl(1) + 40 * scale_y, sprintf('%.1fcm/y', v_legend / 10), 'FontSize', 10);
                else
                    ht = text(xl(1) + 30 * scale_x, yl(1) + 40 * scale_y, sprintf('%.1fcm', v_legend / 10), 'FontSize', 10);
                end
                ht.UserData = struct('keep_color', true);
                ah = arrow(xl(1) + 30 * scale_x, yl(1) + 30 * scale_y, v_legend / lim_vel * max_size_x, 0, h_scale, Core_UI.BLACK);
            end
            ax.Units = old_state;
            max_size = [max_size_x max_size_y];
        end
        
        function [ah_list, ph_list] = insertArrows(fh, id_ref, x, y, v_enu, max_size, lim_vel, cvel_lim, cmap, arrow_type, flag_dot)
            % Insert velocity arrows in a map  of stations
            %
            % SYNTAX
            %   [ah_list, ph_list] = Coordinates.insertArrows(fh, id_ref, x, y, v_enu, vel_lim, cvel_lim, cmap, arrow_type, flag_dot)
            
            if nargin < 11 || isempty(flag_dot)
                flag_dot = false;
            end
            max_size_x = max_size(1);
            max_size_y = max_size(2);
            
            h_scale = 18;
            ax = setAxis(fh);
            
            % Draw the arrows -----------------------------------------
            
            ah_list = cell(length(x),1);
            ph_list = cell(length(x),1);
            for r = 1 : length(x)
                if isempty(id_ref) || (r ~= id_ref)
                    if (strcmp(arrow_type, 'up'))
                        modulus = abs(v_enu(r,3)) * 1e3;
                        id_col = min(size(cmap,1), max(1, ceil(size(cmap,1) * (modulus / cvel_lim(end)))));
                        flag_oversize = (abs(v_enu(r,3) * 1e3 / lim_vel) > 1);
                        if flag_oversize
                            ah = arrow(x(r), y(r), 0,sign(v_enu(r,3)) * max_size_y, h_scale, cmap(id_col, :), flag_oversize);
                        else
                            ah = arrow(x(r), y(r), ...
                                0, ...
                                max(-1,min(1, v_enu(r,3) * 1e3 / lim_vel)) * max_size_y, ...
                                h_scale, cmap(id_col, :));
                        end
                        ah_list{r} = ah;
                    else
                        modulus = sqrt(v_enu(r,1).^2 + v_enu(r,2).^2) * 1e3;
                        id_col = min(size(cmap,1), max(1, ceil(size(cmap,1) * (modulus / cvel_lim(1)))));
                        flag_oversize = (abs(v_enu(r,1) * 1e3 / lim_vel) > 1) || (abs(v_enu(r,2) * 1e3 / lim_vel) > 1);
                        if flag_oversize
                            sf = max(abs(v_enu(r,1) * 1e3 / lim_vel), abs(v_enu(r,2) * 1e3 / lim_vel));
                        else
                            sf = 1;
                        end
                        [ah, ph] = arrow(x(r), y(r), ...
                            (v_enu(r,1) * 1e3 / lim_vel / sf) * max_size_x, ...
                            (v_enu(r,2) * 1e3 / lim_vel / sf) * max_size_y, ...
                            h_scale, cmap(id_col, :), flag_oversize, flag_dot);
                        
                        ah_list{r} = ah;
                        ph_list{r} = ph;
                    end
                end
            end
        end
        
        function mouseMove(~, ~, ax, p_pos, ph)
            % Auxilliary function to showMap
            % Display info-boxes on mouse move
            %
            % SYNTAX
            %   mouseMove(~, ~, ax, p_pos, ph);
            %
            % SEE ALSO
            %   showMap
            
            % get x and y position of cursor
            x_pos = ax.CurrentPoint(1,1);
            y_pos = ax.CurrentPoint(1,2);
            
            old_state = ax.Units;
            ax.Units = 'Points';
            if isempty(ax.UserData)
                ax.UserData = struct('active_labels', []);
            elseif isstruct(ax.UserData)
                if not(isfield(ax.UserData, 'active_labels'))
                    ax.UserData.active_labels = [];
                end
            end
            % Distance between the cursor and the points
            d = sqrt(((p_pos(:,1) - x_pos) / diff(xlim) * ax.Position(3)).^2 + ((p_pos(:,2) - y_pos) / diff(ylim) * ax.Position(4)).^2);
            for l = setdiff(ax.UserData.active_labels, find(d < 30))
                if not(isempty(l))
                    ph(l).UserData.label.Visible = 'off';
                    ax.UserData.active_labels = setdiff(ax.UserData.active_labels, l);
                end
            end
            [d, id_min] = min(d);
            id_ok = (verLessThan('matlab', '99.9') && isfield(ph(id_min).UserData,'label')) || isfield(ph(id_min).UserData,'datatip');
            if id_ok && d < 6
                if verLessThan('matlab', '99.9')
                    ph(id_min).UserData.label.Visible = 'on';
                    scale_x = diff(xlim) / ax.Position(3);
                    txt = ph(id_min).UserData.label;
                    txt.Position(1) = p_pos(id_min,1) - txt.Extent(3) - 20 * scale_x;
                    txt.Position(3) = 1; % on top
                else
                    % data tip are disable untill matlab 3000 they are too slow!
                    ph(id_min).UserData.datatip.Visible = 'on';
                end
            elseif d > 20
                id_min = [];
            end
            % Close all the other info text
            % {id_min '-' ax.UserData.active_labels}
            for l = setdiff(ax.UserData.active_labels, id_min)
                if not(isempty(l))
                    if verLessThan('matlab', '99.9')
                        ph(l).UserData.label.Visible = 'off';
                    else
                        ph(l).UserData.datatip.Visible = 'off';
                    end
                    ax.UserData.active_labels = setdiff(ax.UserData.active_labels, l);
                end
            end
            if id_ok && d < 6
                % Add this new label to the list of open labels
                ax.UserData.active_labels = unique([ax.UserData.active_labels id_min]);
            end
            ax.Units = old_state;
            % display distances in textbox (debug)
            % ax.Title.String = "d: "+d+" id min: "+id_min;
        end
        
        function fixLabels(lbl_pos, lbl_handle)
            % Auxilliary function to showMap
            % Reposition Labels trying not to overlap them
            % Brute force, but somehow working
            %
            % INPUT:
            %   lbl_pos       positions of the labels [X, Y, W, H]
            %                 - X,Y are the coordinates to the lower left corner
            %                 - W,H are width and height
            %   lbl_handle    handle to the lables (text objects) [n x 1]
            %                 (or [n x 2] if handles to labels background are provided)
            %
            % SYNTAX
            %   fixLabel(lbl_pos, lable_handle);
            %
            % SEE ALSO
            %   showMap
            
            % If I have more than one label
            if size(lbl_pos, 1) > 1
                lbl_scale = 5e3/max((diff(minMax(lbl_pos(:, 1))) + lbl_pos(1,3)), (diff(minMax(lbl_pos(:, 2))) + lbl_pos(1,3)));
                x0_lbl = (lbl_pos(:,1) - lbl_pos(:,3)) * lbl_scale;
                y0_lbl = (lbl_pos(:,2) - 0.5 * lbl_pos(:,4)) * lbl_scale;
                x2_lbl = (lbl_pos(:,1) + lbl_pos(:,3)) * lbl_scale;
                y2_lbl = (lbl_pos(:,2) + 1.5 * lbl_pos(:,4)) * lbl_scale;
                b_area = [min(x0_lbl) min(y0_lbl) max(x2_lbl) max(y2_lbl)];
                
                x0_lbl = round(x0_lbl - b_area(1)) +1;
                y0_lbl = round(y0_lbl - b_area(2)) +1;
                
                x2_lbl = round(x2_lbl - b_area(1)) +1;
                y2_lbl = round(y2_lbl - b_area(2)) +1;
                
                w_lbl = (x2_lbl - x0_lbl)/2;
                h_lbl = (y2_lbl - y0_lbl)/2;
                
                b_area = round(b_area - [b_area(1:2) b_area(1:2)]) + 1;
                box = zeros(b_area(4), b_area(3));
                % DEBUG: final_box = zeros(b_area(4), b_area(3));
                for i = 1:size(lbl_handle, 1)
                    box(y0_lbl(i):y2_lbl(i), x0_lbl(i):x2_lbl(i)) = box(y0_lbl(i):y2_lbl(i), x0_lbl(i):x2_lbl(i)) + 1;
                end
                
                if any(box(:) > 1) % I might have overlap of boxes
                    id_min = 1;
                    missing_lbl = 1 : size(lbl_handle, 1);    % List of labels to fit
                    while not(isempty(missing_lbl))
                        score = 0;
                        while score == 0 % repeat untill one lbl have been selected
                            box_score = zeros(numel(missing_lbl),1);
                            for i = 1:numel(missing_lbl) % for each label to fit compute the score (label with higher "free" space)
                                box_score(i) = sum(sum(box(y0_lbl(i) : y2_lbl(i), x0_lbl(i) : x2_lbl(i)) <= id_min));
                            end
                            [score, id_next] = max(box_score);
                            if score == 0
                                id_min = id_min + 1; % I only have overlapped labels
                            end
                        end
                        lbl_id = missing_lbl(id_next);
                        % "Remove" the current box from the conflict box
                        % I only remove 0.5 because I like the idea not to overlap labels in the bounding boxes around points
                        box(y0_lbl(lbl_id) : y2_lbl(lbl_id), x0_lbl(lbl_id):x2_lbl(lbl_id)) = box(y0_lbl(lbl_id) : y2_lbl(lbl_id), x0_lbl(lbl_id):x2_lbl(lbl_id)) - 0.25;
                        % Chose the best position against 8 possible box locations:
                        % 3 8 4
                        % 7 x 5     x is the location of the node to label
                        % 2 6 1
                        pos_limits = round([ ...
                            x0_lbl(lbl_id) + w_lbl(lbl_id)     y0_lbl(lbl_id); ...
                            x0_lbl(lbl_id)                     y0_lbl(lbl_id); ...
                            x0_lbl(lbl_id)                     y0_lbl(lbl_id) + h_lbl(lbl_id); ...
                            x0_lbl(lbl_id) + w_lbl(lbl_id)     y0_lbl(lbl_id) + h_lbl(lbl_id); ...
                            x0_lbl(lbl_id) + w_lbl(lbl_id)     y0_lbl(lbl_id) + h_lbl(lbl_id)/2; ...
                            x0_lbl(lbl_id) + w_lbl(lbl_id)/2   y0_lbl(lbl_id); ...
                            x0_lbl(lbl_id)                     y0_lbl(lbl_id) + h_lbl(lbl_id)/2; ...
                            x0_lbl(lbl_id) + w_lbl(lbl_id)/2   y0_lbl(lbl_id) + h_lbl(lbl_id)]);
                        
                        pos_limits = round([pos_limits(:,1), pos_limits(:,2), min(size(box,2), pos_limits(:,1) + w_lbl(lbl_id)), min(size(box,1), pos_limits(:,2)  + h_lbl(lbl_id))]);
                        i = 1;
                        score = 0;
                        box_score = zeros(8, 1);
                        box_score_sum = zeros(8, 1);
                        favourite_pos = [5 7 6 8 1:4];
                        while i <= 8 && score == 0 % for each label to fit compute the score (label with higher "free" space)
                            p = favourite_pos(i);
                            if all(box(pos_limits(p,2) : pos_limits(p,4), pos_limits(p,1) : pos_limits(p,3)) < 1)
                                % we have our label position!
                                box_score(p) = inf;
                                score = 1;
                            else
                                box_score(p) = sum(sum(box(pos_limits(p,2) : pos_limits(p,4), pos_limits(p,1) : pos_limits(p,3)) <= id_min));
                                box_score_sum(p) = sum(sum(box(pos_limits(p,2) : pos_limits(p,4), pos_limits(p,1) : pos_limits(p,3))));
                            end
                            i = i + 1;
                        end
                        
                        [score] = max(box_score); % We have found our candidate to position
                        id_pos = find(box_score == score); % multiple candidates
                        [score] = min(box_score_sum(id_pos));
                        id_pos = id_pos(find(box_score_sum(id_pos) == score, 1, 'first')); % Get the first good
                        
                        % Save in the overlap box that some areas are now occupied by a label
                        box(pos_limits(id_pos,2) : pos_limits(id_pos,4), pos_limits(id_pos,1) : pos_limits(id_pos,3)) = 1000;
                        % Save in the final map where the labels are located
                        % DEBUG: final_box(pos_limits(id_pos,2) : pos_limits(id_pos,4), pos_limits(id_pos,1) : pos_limits(id_pos,3)) = final_box(pos_limits(id_pos,2) : pos_limits(id_pos,4), pos_limits(id_pos,1) : pos_limits(id_pos,3)) + 1;
                        
                        % 3 8 4
                        % 7 x 5     x is the location of the node to label
                        % 2 6 1
                        w_tmp = lbl_pos(lbl_id,3);
                        h_tmp = lbl_pos(lbl_id,4);
                        x0_tmp = lbl_pos(lbl_id,1) - w_tmp;
                        y0_tmp = lbl_pos(lbl_id,2) - h_tmp/2;
                        h_margin = h_tmp/4;
                        w_margin = h_tmp/4;
                        
                        switch id_pos
                            case 1, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp + w_tmp - w_margin    y0_tmp - h_margin];
                            case 2, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp         + w_margin    y0_tmp - h_margin];
                            case 3, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp         + w_margin    y0_tmp + h_tmp + h_margin];
                            case 4, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp + w_tmp - w_margin    y0_tmp + h_tmp + h_margin];
                            case 5, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp + w_tmp               y0_tmp + h_tmp/2];
                            case 6, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp + w_tmp/2             y0_tmp - h_margin];
                            case 7, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp                       y0_tmp + h_tmp/2];
                            case 8, lbl_handle(lbl_id, 1).Position(1:2) = [x0_tmp + w_tmp/2             y0_tmp + h_tmp + h_margin];
                        end
                        if size(lbl_handle, 2) == 2
                            lbl_handle(lbl_id, 2).Position(1:2) = lbl_handle(lbl_id, 1).Position(1:2);
                        end
                        missing_lbl = setdiff(missing_lbl, lbl_id);
                    end
                end
                
                % figure; imagesc(final_box);
                % figure; imagesc(flipud(box)); caxis([0 5]);
            end
        end
    end
    
    % =========================================================================
    %%    PATHS
    % =========================================================================
    
    methods (Access = 'public')
        
        function out_file_path = getOutPath(this, out_file_prefix)
            % Get the path to ag eneriv coordinate file (noextension)
            %
            % SYNTAX
            %   out_file_path = this.getOutPath(<out_file_prefix>)
            
            state = Core.getState();
            if nargin < 2 || isempty(out_file_prefix)
                out_file_prefix = strrep([state.getPrjName '_'], ' ', '_');
            end
            % Add the folder if not present
            if sum(out_file_prefix == filesep) == 0
                out_dir = state.getOutDir();
                out_file_prefix = fullfile(out_dir, out_file_prefix);
            end
            out_file_path = out_file_prefix;
        end
        
        function out_file_name = getMatOutPath(this, out_file_prefix)
            % Get the path to the coordinate file
            %
            % SYNTAX
            %   out_file_path = this.getCooOutPath(<out_file_prefix>)
            
            % if numel(this) == 1
            %     % add the name to the filename
            %     if (nargin == 2)
            %         out_file_name = strrep([this.getOutPath(out_file_prefix) this.name_v3 '_coo.mat'], ' ', '_');
            %     else
            %         out_file_name = strrep([this.getOutPath() this.name_v3 '_coo.mat'], ' ', '_');
            %     end
            % else
                if (nargin == 2)
                    out_file_name = strrep([this.getOutPath(out_file_prefix) 'coo.mat'], ' ', '_');
                else
                    out_file_name = strrep([this.getOutPath() 'coo.mat'], ' ', '_');
                end
           % end
        end
        
        function out_file_name = getCooOutPath(this, out_file_prefix)
            % Get the path to the coordinate file
            %
            % SYNTAX
            %   out_file_path = this.getCooOutPath(<out_file_prefix>)
            
            if (nargin == 2)
                out_file_name = strrep([this.getOutPath(out_file_prefix) this.name '.coo'], ' ', '_');
            else
                out_file_name = strrep([this.getOutPath() this.name '.coo'], ' ', '_');
            end
        end
        
        function out_file_name = getBernyOutPath(this, out_file_prefix)
            % Get the path to the coordinate file
            %
            % SYNTAX
            %   out_file_path = this.getCooOutPath(<out_file_prefix>)
            
            if (nargin == 2)
                out_file_name = strrep([this.getOutPath(out_file_prefix) this.name '.brn'], ' ', '_');
            else
                out_file_name = strrep([this.getOutPath() this.name '.brn'], ' ', '_');
            end
        end
        
        function out_file_name = getReflectorTxtOutPath(this, out_file_prefix)
            % Get the path to the reflector file
            %
            % SYNTAX
            %   out_file_path = this.getReflectorTxtOutPath(<out_file_prefix>)
            
            if (nargin < 2) || isempty(out_file_prefix)
                out_file_name = strrep([this.getOutPath() this.name '_rfx.txt'], ' ', '_');
            else
                out_file_name = strrep([this.getOutPath(out_file_prefix) this.name '_rfx.txt'], ' ', '_');
            end
        end
        
        function out_file_name = getPrjReflectorTxtOutPath(this, out_file_prefix)
            % Get the path to the reflector file for the mean value of the project
            %
            % SYNTAX
            %   out_file_path = this.getPrjReflectorTxtOutPath(<out_file_prefix>)
            
            if (nargin == 2)
                out_file_name = strrep([this.getOutPath(out_file_prefix) 'merged_rfx.txt'], ' ', '_');
            else
                out_file_name = strrep([this.getOutPath() 'merged_rfx.txt'], ' ', '_');
            end
        end
        
        function out_file_path = exportAsMat(coo_list, out_file_path, flag_rfl_only)
            % Export as coo file (progressive appended file)
            % Any new entry is inserted sorted in the file
            %
            % INPUT
            %   out_file_name   full path of the filename (as default exported into outDir with the name of the coo)
            %
            % SYNTAX
            %   out_file_path = coo.exportAsMat(<out_file_name>)
            
            if nargin < 2 || isempty(out_file_path)
                out_file_path = coo_list.getMatOutPath();
            end
            if nargin < 3 || isempty(flag_rfl_only)
                flag_rfl_only = false;
            end
            log  = Logger.getInstance;
            log.addMarkedMessage(sprintf('Updating coordinates to %s', out_file_path));
            
            [base, fname, fext] = fileparts(out_file_path);
            if ~isempty(base) && not(exist(base, 'dir'))
                try
                    mkdir(base);
                catch ex
                    Core_UI.printEx(ex);
                end
            end
            flag_old_coo = false;
            
            try
                load(out_file_path, 'coo');
                flag_old_coo = not(all(coo.isEmpty));
            catch
                % There are no old coordinates
            end
            
            if flag_old_coo
                % Import old coordinates
                for r = 1 : numel(coo)
                    name_v3 = coo(r).getNameV3;
                    [tmp_coo, id_coo] = coo_list.get(name_v3);
                    if not(isempty(id_coo))
                        coo(r).append(tmp_coo(1), flag_rfl_only);
                        coo_list(id_coo) = [];
                    end
                end
                coo = [coo coo_list];
            else
                coo = coo_list;
            end
            for c = numel(coo):-1:1
                if ~any(coo(c).xyz(:))
                    log.addWarning(sprintf('A coordinate named "%s" does not have valid coordinates, it will not be exported', coo(c).getNameV3));
                    coo(c) = [];
                end
            end
            Coordinates.safeSave(out_file_path, coo);
        end
            
        function exportAsCoo(coo_list, out_file_name, flag_latlonup)
            % Export as coo file (progressive appended file)
            % Any new entry is inserted sorted in the file
            %
            % INPUT
            %   out_file_namArrowse   full path of the filename (as default exported into outDir with the name of the coo)
            %
            % SYNTAX
            %   coo.exportAsCoo(<out_file_name>)
            
            for coo = coo_list(:)'
                now_time = GPS_Time.now();
                if nargin < 2 || isempty(out_file_name)
                    out_file_name = coo.getCooOutPath();
                end
                if nargin < 3 || isempty(flag_latlonup)
                    flag_latlonup = false;
                end
                log  = Logger.getInstance;
                log.addMarkedMessage(sprintf('Updating coordinates to %s', out_file_name));
                try
                    
                    if exist(out_file_name, 'file') == 2
                        % Read and append
                        [txt, lim] = Core_Utils.readTextFile(out_file_name, 3);
                        if isempty(lim)
                            file_ok = false;
                            timestamp = [];
                        else
                            % Verify the file version (it should match 1.0):
                            id_ver = find(txt(lim(:,1) + 1) == 'F'); % +FileVersion
                            version_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), ['(?<=FileVersion[ ]*: )' coo.VERSION], 'once')));
                            if not(version_ok)
                                log = Logger.getInstance;
                                log.addWarning(sprintf('"%s" is in an older format', out_file_name));
                            end
                            file_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), '(?<=FileVersion[ ]*: )1.', 'once')));
                            
                            % Data should be present
                            timestamp = [];
                            if file_ok
                                id_len_ok = find(lim(:,3) >= 9);
                                try
                                    data_start = id_len_ok(find(txt(lim(id_len_ok,1) + 9) == 't') + 1); % +DataStart
                                    id_len_ok = find(lim(:,3) >= 8);
                                    data_stop = id_len_ok(find(txt(lim(id_len_ok,1) + 7) == 'd') -1); % +DataStop
                                    if isempty(data_stop)
                                        data_stop = size(lim, 1);
                                    end
                                    if isempty(data_start)
                                        file_ok = false;
                                    else
                                        id_data = lim(data_start:data_stop,1);
                                        % Read old timestamps
                                        timestamp = datenum(txt(repmat(id_data, 1, 19) + repmat(0:18, numel(id_data), 1)), 'yyyy-mm-dd HH:MM:SS');
                                    end
                                catch
                                    file_ok = false;
                                    timestamp = [];
                                end
                            else
                                data_start = 0;
                            end
                            
                            % Check description field
                            % use the old one for the file
                            if file_ok
                                id_descr = find(txt(lim(:,1) + 4) == 'c'); % + Description
                                if isempty(id_descr)
                                    file_ok = false;
                                else
                                    str_tmp = sprintf('%s\n', txt(lim(id_descr,1):lim(id_descr,2))); % Keep the description of the old file
                                end
                            end
                        end
                    else
                        file_ok = false;
                        timestamp = [];
                    end
                    
                    if not(file_ok)
                        str_tmp = sprintf('+Description    : XYZ Position file generated on %s\n', now_time.toString('dd-mmm-yyyy HH:MM'));
                    end
                    [name, descr] = coo.getName();
                    str_tmp = sprintf('%s+LastChange     : %s\n', str_tmp, now_time.toString('dd-mmm-yyyy HH:MM'));
                    str_tmp = sprintf('%s+Software       : goGPS\n', str_tmp);
                    str_tmp = sprintf('%s+Version        : %s\n', str_tmp, Core.APP_VERSION);
                    str_tmp = sprintf('%s+FileVersion    : %s\n', str_tmp, coo.VERSION);
                    str_tmp = sprintf('%s+MonitoringPoint: %s\n', str_tmp, coo.name);
                    str_tmp = sprintf('%s+NameV3         : %s\n', str_tmp, coo.name_v3);
                    str_tmp = sprintf('%s+LongName       : %s\n', str_tmp, coo.description);
                    str_tmp = sprintf('%s+SensorType     : GNSS\n', str_tmp);
                    str_tmp = sprintf('%s+SensorName     : GNSS\n', str_tmp);
                    if flag_latlonup
                        str_tmp = sprintf('%s+DataScale      : m, deg\n', str_tmp);
                    else
                        str_tmp = sprintf('%s+DataScale      : m\n', str_tmp);
                    end
                    str_tmp = sprintf('%s+DataScale Cov  : mm^2\n', str_tmp);
                    str_tmp = sprintf('%s+DataRate       : %f s\n', str_tmp, coo.getRate);
                    str_tmp = sprintf('%s+DataType       :\n', str_tmp);
                    str_tmp = sprintf('%s -00            : timeStamp\n', str_tmp);
                    str_tmp = sprintf('%s -01            : exportTime\n', str_tmp);
                    str_tmp = sprintf('%s -02            : x\n', str_tmp);
                    str_tmp = sprintf('%s -03            : y\n', str_tmp);
                    str_tmp = sprintf('%s -04            : z\n', str_tmp);
                    str_tmp = sprintf('%s -05            : Cxx\n', str_tmp);
                    str_tmp = sprintf('%s -06            : Cyy\n', str_tmp);
                    str_tmp = sprintf('%s -07            : Czz\n', str_tmp);
                    str_tmp = sprintf('%s -08            : Cxy\n', str_tmp);
                    str_tmp = sprintf('%s -09            : Cxz\n', str_tmp);
                    str_tmp = sprintf('%s -10            : Cyz\n', str_tmp);
                    str_tmp = sprintf('%s -11            : nEpochs\n', str_tmp);
                    str_tmp = sprintf('%s -12            : nObs\n', str_tmp);
                    str_tmp = sprintf('%s -13            : initialSigma0\n', str_tmp);
                    str_tmp = sprintf('%s -14            : sigma0\n', str_tmp);
                    str_tmp = sprintf('%s -15            : fixingRatio\n', str_tmp);
                    str_tmp = sprintf('%s -16            : obsRate\n', str_tmp);
                    str_tmp = sprintf('%s -17            : cooType\n', str_tmp);
                    str_tmp = sprintf('%s -18            : masterName\n', str_tmp);
                    if flag_latlonup
                        str_tmp = sprintf('%s -19            : lat\n', str_tmp);
                        str_tmp = sprintf('%s -20            : lon\n', str_tmp);
                        str_tmp = sprintf('%s -21            : up (ellipsoidal WGS84)\n', str_tmp);
                    end
                    str_tmp = sprintf('%s+DataStart\n', str_tmp);
                    
                    % Append New
                    e = 1; % old epoch
                    [~, id_time] = sort(coo.time.getMatlabTime);
                    if flag_latlonup
                        [lat, lon, up] = coo.getGeodetic();
                    end
                    for i = id_time(:)'
                        cur_time = round(coo.time.getEpoch(i).getMatlabTime*86400)/86400;
                        while e <= numel(timestamp) && (cur_time - 1e-5 > timestamp(e))
                            old_line = txt(lim(data_start + (e-1),1):lim(data_start + (e-1),2));
                            str_tmp = sprintf('%s%s\n', str_tmp, old_line);
                            e = e + 1;
                        end
                        try
                            time = coo.time.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS');
                            xyz = coo.xyz(i,:);
                            if isempty(coo.Cxx) || (i > size(coo.Cxx,3))
                                cov = nan(3,3);
                            else
                                cov = coo.Cxx(:,:,i)*1e6;
                            end
                            try
                                n_epo = coo.info.n_epo(i);
                            catch
                                n_epo = nan;
                            end
                            try
                                n_obs = coo.info.n_obs(i);
                            catch
                                n_obs = nan;
                            end
                            try
                                fix_ratio = coo.info.fixing_ratio(i);
                            catch
                                fix_ratio = nan;
                            end
                            try
                                rate = coo.info.rate(i);
                            catch
                                rate = nan;
                            end
                            try
                                s0_ip = coo.info.s0_ip(i);
                            catch
                                s0_ip = nan;
                            end
                            try
                                s0 = coo.info.s0(i);
                            catch
                                s0 = nan;
                            end
                            try
                                coo_type = char(coo.info.coo_type(i));
                            catch
                                coo_type = 'U';
                            end
                            try
                                master_name = coo.info.master_name(i);
                                [master_name] = Core_Utils.markerCode2MarkerName(master_name); % Convert using the position of the local receiver (I hope they are in the same Country)
                            catch
                                master_name = coo.getNameV3; % Convert marker_v3
                            end
                            if flag_latlonup                                
                                str_tmp = sprintf('%s%s;%s;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%d;%d;%.3f;%.4f;%.2f;%d;%c;%s;%.9f;%.9f;%.4f\n', str_tmp, time, now_time.toString('yyyy-mm-dd HH:MM:SS'), ...
                                xyz(1), xyz(2), xyz(3), ...
                                cov(1,1), cov(2,2), cov(3,3), cov(1,2), cov(1,3), cov(2,3), ...
                                n_epo, ...
                                n_obs, ...
                                s0_ip, ...
                                s0, ...
                                fix_ratio, ...
                                rate, ...
                                char(coo_type), ...
                                master_name, ...
                                lat(i)/pi*180, lon(i)/pi*180, up(i));
                            else
                                str_tmp = sprintf('%s%s;%s;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%.4f;%d;%d;%.3f;%.4f;%.2f;%d;%c;%s\n', str_tmp, time, now_time.toString('yyyy-mm-dd HH:MM:SS'), ...
                                xyz(1), xyz(2), xyz(3), ...
                                cov(1,1), cov(2,2), cov(3,3), cov(1,2), cov(1,3), cov(2,3), ...
                                n_epo, ...
                                n_obs, ...
                                s0_ip, ...
                                s0, ...
                                fix_ratio, ...
                                rate, ...
                                char(coo_type), ...
                                master_name);
                            end
                        catch ex
                            % There is an inconsistency with the entry
                            % could not add this epoch
                            log.addWarning(sprintf('There is a corrupted coordinate in "%s"', name));
                        end
                        % Skip recomputed old epochs
                        while e <= numel(timestamp) && (abs(cur_time - timestamp(e)) < 1e-5)
                            e = e + 1;
                        end
                    end
                    %  Insert old epochs not yet recomputed
                    while e <= numel(timestamp)
                        old_line = txt(lim(data_start + (e-1),1):lim(data_start + (e-1),2));
                        str_tmp = sprintf('%s%s\n', str_tmp, old_line);
                        e = e +1;
                    end
                    fid = fopen(out_file_name, 'Wb');
                    fprintf(fid, str_tmp);
                    fprintf(fid, '+DataEnd\n');
                    fclose(fid);
                    log.addStatusOk(sprintf('Exporting completed successfully'));
                catch ex
                    Core_Utils.printEx(ex);
                    log.addError(sprintf('Exporting failed'));
                end
		out_file_name = '';
            end
        end
        
        function out_file_name = exportAsBerny(this, out_file_name)
            % Export as CRD and OUT file as Bernese does
            %
            % INPUT
            %   out_file_name      full path of the filename (as default exported into outDir with the name of the coo)
            %
            % SYNTAX
            %   coo.exportAsBerny(<out_file_name>)
            %   coo.exportAsBerny(<out_file_name>, new_ref_name, new_fixed_xyz, keep_orphans)
            
            if ~isempty(this) && ~all(this.isEmpty)
                now_time = GPS_Time.now();
                if nargin < 2 || isempty(out_file_name)
                    out_file_name = this.getBernyOutPath();
                end
                %remove extension if present
                [base, fname, fext] = fileparts(out_file_name);
                if not(isempty(base)) && not(exist(base, 'dir'))
                    try
                        mkdir(base);
                    catch ex
                        Core_UI.printEx(ex);
                    end
                end
                out_file_name = fullfile(base, fname);
                clear base name ext;

                log  = Logger.getInstance;
                log.addMarkedMessage(sprintf('Berny Export to "%s"', [out_file_name, '.crd']));
                log.addMarkedMessage(sprintf('Berny Export to "%s"', [out_file_name, '.out']));
                try
                    data_start = 3;
                    if exist([out_file_name, '.crd'], 'file') == 2
                        % Read and append
                        [txt_crd, lim_crd] = Core_Utils.readTextFile([out_file_name, '.crd'], 3);
                        [txt_out, lim_out] = Core_Utils.readTextFile([out_file_name, '.out'], 3);
                        if isempty(lim_crd)
                            file_ok = false;
                            timestamp = [];
                        else
                            % Data should be present
                            file_ok = true;
                            timestamp = [];
                            id_len_ok = find(lim_crd(:,3) >= 9);
                            lim_crd = lim_crd(id_len_ok,:);
                            try
                                % Read old timestamps
                                timestamp = datenum(txt_crd(lim_crd(data_start:end,1) + repmat([16:25 28:36],size(lim_crd,1)-2,1)), 'yyyy-mm-dd HH:MM:SS');
                            catch
                                file_ok = false;
                                timestamp = [];
                            end
                        end
                    else
                        file_ok = false;
                        timestamp = [];
                    end

                    % Write header OUT
                    str_out = '';
                    str_out = sprintf('%swwww-d yyyy-ddd yyyy-mm-dd s   start    end    n epochs    n C1     n C2     n L1     n L2    RMSnet (mm) Rate (s) \n',str_out);
                    str_out = sprintf('%s------+--------+----------+-+--------+--------+---------+--------+--------+--------+--------+------------+--------+\n',str_out);

                    % Write header CRD
                    str_crd = '';
                    str_crd = sprintf('%swwww-d yyyy-ddd yyyy-mm-dd s   start    end    MAST  CRMS        X (m)          Y (m)            Z (m)         EAST (m)       NORTH (m)         UP (m)     F  sE(mm)   sN(mm)   sU(mm)  s3D (mm) \n',str_crd);
                    str_crd = sprintf('%s------+--------+----------+-+--------+--------+----+------+---------------+---------------+---------------+---------------+---------------+---------------+-+--------+--------+--------+--------+\n',str_crd);

                    enu = this.getENU;
                    std_enu = min(this.getStdENU, 99.99*1e-3);
                    std_xyz = min(this.getStdXYZ, 99.99*1e-3);

                    % Get sol_rate
                    sol_rate = this.getRateSolution();

                    % Append New
                    e = 1; % old epoch
                    [~, id_time] = sort(this.time.getMatlabTime);

                    coo_rate = this.getRate;
                    for i = id_time(:)'
                        cur_time = round(this.time.getEpoch(i).getMatlabTime*86400)/86400;
                        while e <= numel(timestamp) && ((cur_time - 1e-5) > (timestamp(e) + (sol_rate / 86400)/2))
                            old_line = txt_crd(lim_crd(data_start + (e-1),1):lim_crd(data_start + (e-1),2));
                            str_crd = sprintf('%s%s\n', str_crd, old_line);
                            old_line = txt_out(lim_out(data_start + (e-1),1):lim_out(data_start + (e-1),2));
                            str_out = sprintf('%s%s\n', str_out, old_line);
                            e = e +1;
                        end
                        try
                            tmp_time = this.time.getEpoch(i);

                            if coo_rate == 86400
                                sss_char = '0';
                            else
                                time_from_start = tmp_time.getMatlabTime;
                                time_from_start = ((time_from_start - floor(time_from_start)) * 86400); % time from the beginning of the day [s]
                                CHAR_LIST = 'abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz01234567';
                                sss_char = CHAR_LIST(floor(time_from_start/coo_rate) + 1);
                            end

                            [year, doy] = tmp_time.getDOY;
                            [week, sow, dow] = tmp_time.getGpsWeek;
                            date_str = tmp_time.toString('yyyy-mm-dd');
                            tmp_time.addIntSeconds(-sol_rate/2);
                            time_start = tmp_time.toString('HH:MM:SS');
                            tmp_time.addIntSeconds(+sol_rate-1);
                            time_stop = tmp_time.toString('HH:MM:SS');
                            xyz = this.xyz(i,:);
                            %if isempty(this.Cxx)
                            %    cov = zeros(3,3);
                            %else
                            %    cov = min(this.Cxx(:,:,i)*1e6, 99.99*1e-3);
                            %end
                            try
                                n_epo = this.info.n_epo(i);
                            catch
                                n_epo = nan;
                            end
                            try
                                n_obs = this.info.n_obs(i);
                            catch
                                n_obs = nan;
                            end
                            try
                                fix_ratio = this.info.fixing_ratio(i);
                            catch
                                fix_ratio = nan;
                            end
                            try
                                rate = this.info.rate(i);
                            catch
                                rate = nan;
                            end
                            try
                                s0_ip = min(this.info.s0_ip(i), 999);
                            catch
                                s0_ip = nan;
                            end
                            if s0_ip > 999
                                s0_ip = 999;
                            end
                            try
                                s0 = min(this.info.s0(i), 999);
                            catch
                                s0 = nan;
                            end
                            if s0 > 1
                                s0 = 0.999;
                            end
                            try
                                coo_type = char(this.info.coo_type(i));
                                if coo_type == 'F'
                                    std_enu(i,:) = nan;
                                    std_xyz(i,:) = nan;
                                end
                            catch
                                coo_type = 'U';
                            end

                            try
                                master_name = this.info.master_name(i);
                                [master_name] = Core_Utils.getMarkerV3(master_name, this); % Convert using the position of the local receiver (I hope they are in the same Country)
                            catch
                                [master_name] = Core_Utils.getMarkerV3(this.name_v3); % Convert marker_v3
                            end
                            master_name = master_name(1:4);

                            str_crd = sprintf('%s%04d-%1d %04d-%03d %s %c %s %s %4s %6.2f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %1c %8.2f %8.2f %8.2f %8.2f\n', str_crd, ...
                                week,dow, ...
                                year,doy, ...
                                date_str, ...
                                sss_char, ...
                                time_start, time_stop, ...
                                master_name, ...
                                s0_ip, ...
                                xyz(1), xyz(2), xyz(3), ...
                                enu(i,1), enu(i,2), enu(i,3), ...
                                char(coo_type), ...
                                std_enu(i,1).*1e3, std_enu(i,2).*1e3, std_enu(i,3).*1e3, ...
                                sqrt(sum(std_xyz(i,:).^2)).*1e3);
                            str_out = sprintf('%s%04d-%1d %04d-%03d %s %c %s %s %8d  %7d  %7d  %7d  %7d %12.2f %8d\n', str_out, ...
                                week,dow, ...
                                year,doy, ...
                                date_str, ...
                                sss_char, ...
                                time_start, time_stop, ...
                                n_epo, 0, 0, n_obs, 0, ...
                                s0 * 1e3, ...
                                round(rate,2));
                        catch ex
                            % There is an inconsistency with the entry
                            % could not add this epoch
                            log.addWarning('There is a corrupted coordinate');
                        end
                        % Skip recomputed old epochs
                        while e <= numel(timestamp) && (abs(cur_time - (timestamp(e) + (sol_rate / 86400)/2)) < 1e-5)
                            e = e +1;
                        end
                    end

                    %  Insert old epochs not yet recomputed
                    while e <= numel(timestamp)
                        old_line = txt_crd(lim_crd(data_start + (e-1),1):lim_crd(data_start + (e-1),2));
                        str_crd = sprintf('%s%s\n', str_crd, old_line);
                        old_line = txt_out(lim_out(data_start + (e-1),1):lim_out(data_start + (e-1),2));
                        str_out = sprintf('%s%s\n', str_out, old_line);
                        e = e +1;
                    end

                    fid_out = fopen([out_file_name '.out'], 'Wb');
                    fprintf(fid_out, str_out);
                    fclose(fid_out);

                    fid_crd = fopen([out_file_name '.crd'], 'Wb');
                    fprintf(fid_crd, str_crd);
                    fclose(fid_crd);

                    log.addStatusOk(sprintf('Exporting completed successfully'));
                catch ex
                    Core_Utils.printEx(ex);
                    log.addError(sprintf('Exporting failed'));
                end
            end
        end
        
        function file_list_path = exportReflectometryTxtInterp(coo_list, out_dir, first_day, last_day, rt_par)
            
            state = Core.getState;
            if nargin < 2 || isempty(out_dir)
                out_dir = state.getOutDir();
                out_dir = fullfile(out_dir, 'refl');
            end
            if nargin < 3 || isempty(first_day)
                first_day = state.getSessionsStart;
            end
            if nargin < 4 || isempty(last_day)
                last_day = state.getSessionsStop;
            end

            file_list_path = {};

            log = Logger.getInstance;
            log.addMarkedMessage(sprintf('Exporting reflector heights'));

            % get a copy, we are going to trim the data
            coo = coo_list.getCopy;
            lim_start = first_day.getCopy; lim_start.addIntSeconds(-2*86400); % trim in window of 10 days
            lim_stop = last_day.getCopy; lim_stop.addIntSeconds(+2*86400); % trim in window of 10 days
            coo.keep(lim_start, lim_stop); % cut to avoid too complex computations
            
            [rh, time_data, err] = coo.getReflectorHeightInterp(86400/4, 900);
            [~, ~, ~, h_masl] = coo.getGeodetic;

            % Get station position, and interpolate
            time_coo = coo.getTime;
            if ~any(rh)
                log.addError('No reflectometry data to export');
                return
            end
            if time_coo.length > 1
                value = interp1(time_coo.getMatlabTime, h_masl, time_data.getMatlabTime, 'spline', 'extrap') - rh;
            else
                value = h_masl - coo.reflectometry.value;
            end

            if time_coo.length == 1
                xyz = coo.xyz;
            else
                x = interp1(time_coo.getMatlabTime, coo.xyz(:,1), time_data.getMatlabTime, 'spline');
                y = interp1(time_coo.getMatlabTime, coo.xyz(:,2), time_data.getMatlabTime, 'spline');
                z = interp1(time_coo.getMatlabTime, coo.xyz(:,3), time_data.getMatlabTime, 'spline');
                xyz = [x y z];
            end

            %% Export a file for each day files will be ovewritten with the newest interpolation
            first_day = floor(first_day.getMatlabTime);
            last_day = floor(last_day.getMatlabTime);
            time_str = time_data.toString('yyyy-mm-dd HH:MM:SS');
            time_data = time_data.getMatlabTime;
            for d = first_day:last_day
                e_start  = find(time_data >= d, 1, 'first');
                e_stop  = find(time_data < d+1, 1, 'last');
                day = GPS_Time(d);
                [year, doy] = day.getDOY;
                cur_out_dir = fullfile(out_dir, num2str(year));
                if ~exist(cur_out_dir,'dir')
                    mkdir(cur_out_dir);
                end
                out_file_path = fullfile(cur_out_dir, [coo.getNameV3,'_', num2str(year), num2str(doy,'%03d'), '_rt_rfx.txt']);

                % Start exporting operation:
                % Prepare header
                str_tmp = sprintf('+Description    : Approximate reflector height [m] on the geoid WGS84 generated on %s\n', GPS_Time.now.toString('dd-mmm-yyyy HH:MM'));

                [name, descr] = coo.getName();
                str_tmp = sprintf('%s+LastChange     : %s\n', str_tmp, GPS_Time.now.toString('dd-mmm-yyyy HH:MM'));
                str_tmp = sprintf('%s+Software       : Breva\n', str_tmp);
                str_tmp = sprintf('%s+Version        : %s\n', str_tmp, Core.APP_VERSION);
                str_tmp = sprintf('%s+FileVersion    : %s\n', str_tmp, coo.VERSION);
                str_tmp = sprintf('%s+MonitoringPoint: %s\n', str_tmp, coo.name);
                str_tmp = sprintf('%s+NameV3         : %s\n', str_tmp, coo.name_v3);
                str_tmp = sprintf('%s+LongName       : %s\n', str_tmp, coo.description);
                str_tmp = sprintf('%s+SensorType     : GNSS Reflectometry\n', str_tmp);
                str_tmp = sprintf('%s+SensorName     : GNSS Reflectometry\n', str_tmp);
                str_tmp = sprintf('%s+DataScale      : m (masl)\n', str_tmp);
                str_tmp = sprintf('%s+DataScale std  : m\n', str_tmp);
                str_tmp = sprintf('%s+DataRate       : %f s\n', str_tmp, coo.getRate);
                str_tmp = sprintf('%s+DataType       :\n', str_tmp);
                str_tmp = sprintf('%s -00            : timeStamp\n', str_tmp);
                str_tmp = sprintf('%s -01            : orthometricHeight\n', str_tmp);
                str_tmp = sprintf('%s -02            : x\n', str_tmp);
                str_tmp = sprintf('%s -03            : y\n', str_tmp);
                str_tmp = sprintf('%s -04            : z\n', str_tmp);
                str_tmp = sprintf('%s -05            : std\n', str_tmp);

                str_tmp = sprintf('%s+DataStart\n', str_tmp);

                % prepare data export:
                for e = e_start:e_stop
                    if ~isnan(value(e))
                    str_tmp = sprintf('%s%s;%8.3f;%.4f;%.4f;%.4f;%.4f\n', str_tmp, time_str(e,:), ...
                        value(e), ...
                        xyz(e,1), xyz(e,2), xyz(e,3), ...
                        err(e));
                    end
                end

                file_list_path = [file_list_path; {out_file_path}];
                fid = fopen(out_file_path, 'Wb');
                    fprintf(fid, str_tmp);
                    fprintf(fid, '+DataEnd\n');
                    fclose(fid);
                    log.addStatusOk(sprintf('Exporting of "%s" completed successfully', out_file_path));
            end
        end

        function exportReflectometryTxt(coo_list, out_file_path)
            % Export as Reflector height as Text File
            % Compatibility layer with GeoGuard
            %
            % INPUT
            %   out_file_name      full path of the filename (as default exported into outDir with the name of the coo)
            %
            % SYNTAX
            %   coo.exportReflectometryTxt(<out_file_name>)
            
            for coo = coo_list(:)'
                now_time = GPS_Time.now();
                if nargin < 2 || isempty(out_file_path)
                    out_file_path = coo.getReflectorTxtOutPath();
                end
                log  = Logger.getInstance;
                log.addMarkedMessage(sprintf('Updating reflector height to %s', out_file_path));
                try
                    [file_dir, file_name, file_ext] = fileparts(out_file_path);
                    if exist(file_dir, 'dir') ~= 7
                        mkdir(file_dir);
                    end
                    if exist(out_file_path, 'file') == 2
                        % Read and append
                        [txt, lim] = Core_Utils.readTextFile(out_file_path, 3);
                        if isempty(lim)
                            file_ok = false;
                            timestamp = [];
                        else
                            % Verify the file version (it should match 1.0):
                            id_ver = find(txt(lim(:,1) + 1) == 'F'); % +FileVersion
                            version_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), ['(?<=FileVersion[ ]*: )' coo.VERSION], 'once')));
                            if not(version_ok)
                                log = Logger.getInstance;
                                log.addWarning(sprintf('"%s" is in an older format', out_file_path));
                            end
                            file_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), '(?<=FileVersion[ ]*: )1.', 'once')));
                            
                            % Data should be present
                            timestamp = [];
                            if file_ok
                                id_len_ok = find(lim(:,3) >= 9);
                                try
                                    data_start = id_len_ok(find(txt(lim(id_len_ok,1) + 9) == 't') + 1); % +DataStart
                                    id_len_ok = find(lim(:,3) >= 8);
                                    data_stop = id_len_ok(find(txt(lim(id_len_ok,1) + 7) == 'd') -1); % +DataStop
                                    if isempty(data_stop)
                                        data_stop = size(lim, 1);
                                    end
                                    if isempty(data_start)
                                        file_ok = false;
                                    else
                                        id_data = lim(data_start:data_stop,1);
                                        % Read old timestamps
                                        timestamp = datenum(txt(repmat(id_data, 1, 19) + repmat(0:18, numel(id_data), 1)), 'yyyy-mm-dd HH:MM:SS');
                                    end
                                catch
                                    file_ok = false;
                                    timestamp = [];
                                end
                            else
                                data_start = 0;
                            end
                            
                            % Check description field
                            % use the old one for the file
                            if file_ok
                                id_descr = find(txt(lim(:,1) + 4) == 'c'); % + Description
                                if isempty(id_descr)
                                    file_ok = false;
                                else
                                    str_tmp = sprintf('%s\n', txt(lim(id_descr,1):lim(id_descr,2))); % Keep the description of the old file
                                end
                            end
                        end
                    else
                        file_ok = false;
                        timestamp = [];
                    end
                    
                    if not(file_ok)
                        str_tmp = sprintf('+Description    : Approximate reflector height [m] on the geoid WGS84 generated on %s\n', now_time.toString('dd-mmm-yyyy HH:MM'));
                    end
                    [name, descr] = coo.getName();
                    str_tmp = sprintf('%s+LastChange     : %s\n', str_tmp, now_time.toString('dd-mmm-yyyy HH:MM'));
                    str_tmp = sprintf('%s+Software       : Breva\n', str_tmp);
                    str_tmp = sprintf('%s+Version        : %s\n', str_tmp, Core.APP_VERSION);
                    str_tmp = sprintf('%s+FileVersion    : %s\n', str_tmp, coo.VERSION);
                    str_tmp = sprintf('%s+MonitoringPoint: %s\n', str_tmp, coo.name);
                    str_tmp = sprintf('%s+NameV3         : %s\n', str_tmp, coo.name_v3);
                    str_tmp = sprintf('%s+LongName       : %s\n', str_tmp, coo.description);
                    str_tmp = sprintf('%s+SensorType     : GNSS Reflectometry\n', str_tmp);
                    str_tmp = sprintf('%s+SensorName     : GNSS Reflectometry\n', str_tmp);
                    str_tmp = sprintf('%s+DataScale      : m (masl)\n', str_tmp);
                    str_tmp = sprintf('%s+DataScale std  : m\n', str_tmp);
                    str_tmp = sprintf('%s+DataRate       : %f s\n', str_tmp, coo.getRate);
                    str_tmp = sprintf('%s+DataType       :\n', str_tmp);
                    str_tmp = sprintf('%s -00            : timeStamp\n', str_tmp);
                    str_tmp = sprintf('%s -01            : exportTime\n', str_tmp);
                    str_tmp = sprintf('%s -02            : orthometricHeight\n', str_tmp);
                    str_tmp = sprintf('%s -03            : x\n', str_tmp);
                    str_tmp = sprintf('%s -04            : y\n', str_tmp);
                    str_tmp = sprintf('%s -05            : z\n', str_tmp);
                    str_tmp = sprintf('%s -06            : nObs\n', str_tmp);
                    str_tmp = sprintf('%s -07            : std\n', str_tmp);

                    str_tmp = sprintf('%s+DataStart\n', str_tmp);
                    
                    % Append New
                    e = 1; % old epoch
                    time_refl = coo.getReflectometryTime;
                    [~, id_time] = sort(time_refl.getMatlabTime);
                    for i = id_time(:)'
                        cur_time = round(time_refl.getEpoch(i).getMatlabTime*86400)/86400;
                        while e <= numel(timestamp) && (cur_time - 1e-5 > timestamp(e))
                            old_line = txt(lim(data_start + (e-1),1):lim(data_start + (e-1),2));
                            str_tmp = sprintf('%s%s\n', str_tmp, old_line);
                            e = e + 1;
                        end
                        try
                            [~, ~, ~, h_masl] = coo.getGeodetic;

                            time_coo = coo.getTime;
                            time_data = time_refl.getEpoch(i);
                            if time_coo.length > 1
                                value = interp1(time_coo.getMatlabTime, h_masl, time_data.getMatlabTime, 'spline') - coo.reflectometry.value(i);
                            else
                                value = h_masl - coo.reflectometry.value(i);
                            end

                            if time_coo.length == 1
                                xyz = coo.xyz;
                            else
                                x = interp1(time_coo.getMatlabTime, coo.xyz(:,1), time_data.getMatlabTime, 'spline');
                                y = interp1(time_coo.getMatlabTime, coo.xyz(:,2), time_data.getMatlabTime, 'spline');
                                z = interp1(time_coo.getMatlabTime, coo.xyz(:,3), time_data.getMatlabTime, 'spline');
                                xyz = [x y z];
                            end
                            
                            n_obs = coo.reflectometry.n_obs(i);
                            std_r = coo.reflectometry.std(i);
                            if ~isnan(value)
                                time = time_data.toString('yyyy-mm-dd HH:MM:SS');
                                str_tmp = sprintf('%s%s;%s;%8.3f;%.4f;%.4f;%.4f;%d;%.4f\n', str_tmp, time, now_time.toString('yyyy-mm-dd HH:MM:SS'), ...
                                    value(1), ...
                                    xyz(1), xyz(2), xyz(3), ...
                                    n_obs, ...
                                    std_r);
                            end
                        catch ex
                            % There is an inconsistency with the entry
                            % could not add this epoch
                            log.addWarning(sprintf('There is a corrupted reflector height in "%s"', name));
                        end
                        % Skip recomputed old epochs
                        while e <= numel(timestamp) && (abs(cur_time - timestamp(e)) < 1e-5)
                            e = e + 1;
                        end
                    end
                    %  Insert old epochs not yet recomputed
                    while e <= numel(timestamp)
                        old_line = txt(lim(data_start + (e-1),1):lim(data_start + (e-1),2));
                        str_tmp = sprintf('%s%s\n', str_tmp, old_line);
                        e = e +1;
                    end
                    fid = fopen(out_file_path, 'Wb');
                    fprintf(fid, str_tmp);
                    fprintf(fid, '+DataEnd\n');
                    fclose(fid);
                    log.addStatusOk(sprintf('Exporting completed successfully'));
                catch ex
                    Core_Utils.printEx(ex);
                    log.addError(sprintf('Exporting failed'));
                end
            end
        end
        
        function exportMergedReflectometryTxt(coo_list, out_file_name)
            % Export as Reflector height as Text File
            % Compatibility layer with GeoGuard
            %
            % INPUT
            %   out_file_name      full path of the filename (as default exported into outDir with the name of the coo)
            %
            % SYNTAX
            %   coo.exportReflectometryTxt(<out_file_name>)
            
            [rh_data, time_data, data_sync, n_obs] = coo_list.getMeanReflectionHeight();
            coo = coo_list(1);
            now_time = GPS_Time.now();
            if nargin < 2 || isempty(out_file_name)
                out_file_name = coo(1).getPrjReflectorTxtOutPath();
            end
            log  = Logger.getInstance;
            log.addMarkedMessage(sprintf('Updating reflector height to %s', out_file_name));
            try
                [file_dir, file_name, file_ext] = fileparts(out_file_name);
                if exist(file_dir, 'dir') ~= 7
                    mkdir(file_dir);
                end
                if exist(out_file_name, 'file') == 2
                    % Read and append
                    [txt, lim] = Core_Utils.readTextFile(out_file_name, 3);
                    if isempty(lim)
                        file_ok = false;
                        timestamp = [];
                    else
                        % Verify the file version (it should match 1.0):
                        id_ver = find(txt(lim(:,1) + 1) == 'F'); % +FileVersion
                        version_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), ['(?<=FileVersion[ ]*: )' coo.VERSION], 'once')));
                        if not(version_ok)
                            log = Logger.getInstance;
                            log.addWarning(sprintf('"%s" is in an older format', out_file_name));
                        end
                        file_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), '(?<=FileVersion[ ]*: )1.', 'once')));
                        
                        % Data should be present
                        timestamp = [];
                        if file_ok
                            id_len_ok = find(lim(:,3) >= 9);
                            try
                                data_start = id_len_ok(find(txt(lim(id_len_ok,1) + 9) == 't') + 1); % +DataStart
                                id_len_ok = find(lim(:,3) >= 8);
                                data_stop = id_len_ok(find(txt(lim(id_len_ok,1) + 7) == 'd') -1); % +DataStop
                                if isempty(data_stop)
                                    data_stop = size(lim, 1);
                                end
                                if isempty(data_start)
                                    file_ok = false;
                                else
                                    id_data = lim(data_start:data_stop,1);
                                    % Read old timestamps
                                    timestamp = datenum(txt(repmat(id_data, 1, 19) + repmat(0:18, numel(id_data), 1)), 'yyyy-mm-dd HH:MM:SS');
                                end
                            catch
                                file_ok = false;
                                timestamp = [];
                            end
                        else
                            data_start = 0;
                        end
                        
                        % Check description field
                        % use the old one for the file
                        if file_ok
                            id_descr = find(txt(lim(:,1) + 4) == 'c'); % + Description
                            if isempty(id_descr)
                                file_ok = false;
                            else
                                str_tmp = sprintf('%s\n', txt(lim(id_descr,1):lim(id_descr,2))); % Keep the description of the old file
                            end
                        end
                    end
                else
                    file_ok = false;
                    timestamp = [];
                end
                
                if not(file_ok)
                    str_tmp = sprintf('+Description    : Approximate reflector height [m] on the ellipsoid generated on %s\n', now_time.toString('dd-mmm-yyyy HH:MM'));
                end
                name_list = '';
                namev3_list = '';
                coo_m = coo_list.getMedianPos;
                xyz_med = nan(numel(coo_m),3); % used to compute the baricentric position
                for i = 1:numel(coo_list)
                    name_list = sprintf('%s, %s', name_list, coo_list(i).getName());
                    namev3_list = sprintf('%s, %s', namev3_list, coo_list(i).name_v3);
                    xyz_med(i, :) = coo_m(i).getXYZ;
                end
                if length(name_list) > 3
                    name_list = name_list(3:end);
                    namev3_list = namev3_list(3:end);
                end
                xyz_med = mean(xyz_med,1,'omitnan');
                str_tmp = sprintf('%s+LastChange     : %s\n', str_tmp, now_time.toString('dd-mmm-yyyy HH:MM'));
                str_tmp = sprintf('%s+Software       : Breva\n', str_tmp);
                str_tmp = sprintf('%s+Version        : %s\n', str_tmp, Core.APP_VERSION);
                str_tmp = sprintf('%s+FileVersion    : %s\n', str_tmp, coo.VERSION);
                str_tmp = sprintf('%s+MonitoringPoint: %s\n', str_tmp, name_list);
                str_tmp = sprintf('%s+NameV3         : %s\n', str_tmp, namev3_list);
                str_tmp = sprintf('%s+SensorType     : GNSS Reflectometry\n', str_tmp);
                str_tmp = sprintf('%s+SensorName     : GNSS Reflectometry\n', str_tmp);
                str_tmp = sprintf('%s+DataScale      : m\n', str_tmp);
                str_tmp = sprintf('%s+DataScale std  : m\n', str_tmp);
                str_tmp = sprintf('%s+DataRate       : %f s\n', str_tmp, time_data.getRate);
                str_tmp = sprintf('%s+DataType       :\n', str_tmp);
                str_tmp = sprintf('%s -00            : timeStamp\n', str_tmp);
                str_tmp = sprintf('%s -01            : exportTime\n', str_tmp);
                str_tmp = sprintf('%s -02            : orthometricHeight\n', str_tmp);
                str_tmp = sprintf('%s -03            : x\n', str_tmp);
                str_tmp = sprintf('%s -04            : y\n', str_tmp);
                str_tmp = sprintf('%s -05            : z\n', str_tmp);
                str_tmp = sprintf('%s -06            : nObs\n', str_tmp);
                
                str_tmp = sprintf('%s+DataStart\n', str_tmp);
                
                % Append New
                e = 1; % old epoch
                [~, id_time] = sort(time_data.getMatlabTime);
                for i = id_time(:)'
                    cur_time = round(time_data.getEpoch(i).getMatlabTime*86400)/86400;
                    while e <= numel(timestamp) && (cur_time - 1e-5 > timestamp(e))
                        old_line = txt(lim(data_start + (e-1),1):lim(data_start + (e-1),2));
                        str_tmp = sprintf('%s%s\n', str_tmp, old_line);
                        e = e + 1;
                    end
                    try
                        time = time_data.getEpoch(i).toString('yyyy-mm-dd HH:MM:SS');
                        tot_n_obs = sum(noNaN(n_obs(i,:)));
                        if ~isnan(rh_data(i))
                            str_tmp = sprintf('%s%s;%s;%8.3f;%.4f;%.4f;%.4f;%d\n', str_tmp, time, now_time.toString('yyyy-mm-dd HH:MM:SS'), ...
                                rh_data(i), ...
                                xyz_med(1), xyz_med(2), xyz_med(3), ...
                                tot_n_obs);
                        end
                    catch ex
                        % There is an inconsistency with the entry
                        % could not add this epoch
                        log.addWarning(sprintf('There is a corrupted reflector height in "%s"', name_list));
                    end
                    % Skip recomputed old epochs
                    while e <= numel(timestamp) && (abs(cur_time - timestamp(e)) < 1e-5)
                        e = e + 1;
                    end
                end
                %  Insert old epochs not yet recomputed
                while e <= numel(timestamp)
                    old_line = txt(lim(data_start + (e-1),1):lim(data_start + (e-1),2));
                    str_tmp = sprintf('%s%s\n', str_tmp, old_line);
                    e = e +1;
                end
                fid = fopen(out_file_name, 'Wb');
                fprintf(fid, str_tmp);
                fprintf(fid, '+DataEnd\n');
                fclose(fid);
                log.addStatusOk(sprintf('Exporting completed successfully'));
            catch ex
                Core_Utils.printEx(ex);
                log.addError(sprintf('Exporting failed'));
            end
        end
    end
    
    % =========================================================================
    %%    OPERATIONS
    % =========================================================================
    
    methods (Access = 'public')
        function this = importCoo(this, file_name)
            this = Coordinates.fromCooFile(file_name);
        end
        
        function res = eq(coo1, coo2)
            %%% DESCRIPTION: check if two coordinates are equal
            d = sqrt(sum((coo1.xyz - coo2.xyz).^2, 2));
            res = d < coo1.precision;
        end
        
        function coo_diff = minus(coo1, coo2)
            coo_diff = struct('time', [], 'enu_diff', [], 'xyz_diff', [], 'xyz_model_diff', []);
            if coo1.time.isEmpty || coo2.time.isEmpty
                coo_diff.time = GPS_Time;
                
                id_ok1 = 1 : min(size(coo1.xyz, 1), size(coo2.xyz, 1));
                id_ok2 = id_ok1;
            else
                [common_time, id_ok1, id_ok2] = intersect(coo1.time.getNominalTime.getMatlabTime, coo2.time.getNominalTime.getMatlabTime);
                coo_diff.time = GPS_Time(common_time);
            end
            
            coo_diff.enu_diff = coo1.getElement(id_ok1).getENU - coo2.getElement(id_ok2).getENU;
            coo_diff.xyz_diff = coo1.getElement(id_ok1).xyz - coo2.getElement(id_ok2).xyz;
            try
                coo_diff.xyz_model_diff =  coo1.getElement(id_ok1).xyz_model - coo2.getElement(id_ok2).xyz_model;
            catch ex
                coo_diff.xyz_model_diff = [];
            end
        end
        
        function setNewRef(coo_list, new_ref_name, new_fixed_xyz, keep_orphans)
            % Fix a coordinate to a new value
            %
            % INPUT
            %   new_ref_name       name of the reference coordinate
            %   new_fixed_xyz      new coordinate of the reference
            %   keep_orphans       keep the epoch with master different from the new reference (default)
            %
            % SYNTAX
            %   coo_list.setNewRef(new_ref_name, new_fixed_xyz, keep_orphans);
            if isnumeric(new_ref_name)
                new_ref_id = new_ref_name;
                ref_found = true;
            else
                % find the new reference station
                c = 0;
                ref_found = false;
                while (c < numel(coo_list)) && ~ref_found
                    c = c + 1;
                    tmp = coo_list(c).getNameV3;
                    n_char = min(numel(tmp), numel(new_ref_name)); 
                    if strcmp(new_ref_name(1:n_char), tmp(1:n_char))
                        ref_found = true;
                    end
                end
                if ~ref_found
                    Logger.getInstance.addError('New reference marker not found! Changing reference is not possible.')
                else
                    new_ref_id = c;
                end
            end
            
            if ref_found
                if nargin < 3 || isempty(new_fixed_xyz)
                    try
                        rf = Core.getReferenceFrame;
                        new_fixed_xyz = rf.getCoo(coo_list(new_ref_id).getName, coo_list(new_ref_id).time.last); % fix to the last coordinate in RF
                    catch
                        % any problem with the RF is managed by using median coordinates
                        new_fixed_xyz = [];
                    end
                    if isempty(new_fixed_xyz)
                        % if empty fix to the median value
                        new_fixed_xyz = coo_list(new_ref_id).getMedianPos.getXYZ;
                    end
                end
                
                [~, new_ref_name] = coo_list(new_ref_id).getNameV3;
                
                coo_rate = round(coo_list(new_ref_id).getRate, 3);
                if isempty(coo_rate) || isnan(zero2nan(coo_rate))
                    % try to retrieve rate from the date
                    coo_rate = round(median(diff(coo_list(new_ref_id).time.getRefTime), 'omitnan'), 3);
                end
                if isnan(coo_rate)
                    coo_rate = 1;
                end
                time_ref = coo_list(new_ref_id).time.getRoundedTime(coo_rate);
                time0 = time_ref.first.getMatlabTime;
                tid_ref = time_ref.getRefTime(time0) / coo_rate;
                if any(time0) || (time_ref.length == 0)
                    if (time_ref.length == 0)
                        tid_ref = 1;
                        xyz_corr = round(repmat(new_fixed_xyz, size(coo_list(new_ref_id).xyz,1), 1) - coo_list(new_ref_id).xyz, 6);
                    else
                        xyz_corr = round(repmat(new_fixed_xyz, numel(tid_ref), 1) - coo_list(new_ref_id).xyz, 6);
                    end
                    if not(isempty(coo_list(new_ref_id).xyz_model)) && any(coo_list(new_ref_id).xyz_model(:))
                        xyz_corr_model = round(repmat(new_fixed_xyz, numel(tid_ref), 1) - coo_list(new_ref_id).xyz_model, 6);
                    else
                        xyz_corr_model = 0;
                    end
                    
                    % for each non reference coordinate
                    for c = setdiff(1 : numel(coo_list), new_ref_id)
                        if (time_ref.length == 0)
                            idr = 1;
                        else
                            tid_coo = coo_list(c).time.getRoundedTime(coo_rate).getRefTime(time0)/coo_rate;                        
                            [~, idc, idr] = intersect(tid_coo, tid_ref);
                        end

                        if numel(idr) == 1
                            idc = (1:size(coo_list(c).xyz,1))';
                        end
                        if any(idc)
                            % apply translation
                            coo_list(c).xyz(idc, :) = coo_list(c).xyz(idc, :) + xyz_corr(idr, :);
                            if numel(coo_list(c).xyz) == numel(coo_list(c).xyz_model) && any(xyz_corr_model(:))
                                coo_list(c).xyz_model(idc, :) = coo_list(c).xyz_model(idc, :) + xyz_corr_model(idr, :);
                            end
                            % Covariance propagation with missing cross covariance term
                            
                            if max(idr) <= size(coo_list(new_ref_id).Cxx,3)
                                vcv_ref = coo_list(new_ref_id).Cxx(:, :, idr);
                            else
                                vcv_ref = [];
                            end
                            if isempty(vcv_ref)
                                vcv_ref = zeros(3);
                            end
                            if max(idc) <= size(coo_list(c).Cxx, 3)
                                vcv = coo_list(c).Cxx(:, :, idc);
                                if isempty(vcv)
                                    vcv = nan(3);
                                end
                                if isempty(coo_list(c).Cxx)
                                    coo_list(c).Cxx = nan(3, 3, idc(end));
                                end
                                coo_list(c).Cxx(:, :, idc) = vcv + vcv_ref;
                            end
                            coo_list(c).info.master_name(idc) = new_ref_name;
                            coo_list(c).info.coo_type(idc) = 'G';
                            
                            % remove epochs with no master
                            if nargin > 3 && not(keep_orphans)
                                id_ko = setdiff((1:coo_list(c).time.length)', idc);
                                coo_list(c).rem(id_ko);
                            end
                        end
                    end
                    
                    % Now fix the new reference
                    coo_list(new_ref_id).xyz = repmat(new_fixed_xyz, numel(tid_ref), 1);
                    coo_list(new_ref_id).xyz_model = repmat(new_fixed_xyz, numel(tid_ref), 1);
                    coo_list(new_ref_id).info.master_name = new_ref_name * ones(size(coo_list(new_ref_id).xyz,1), 1,'uint64');
                    coo_list(new_ref_id).info.coo_type = 'F' * char(ones(1, size(coo_list(new_ref_id).xyz,1), 'uint8'));
                    coo_list(new_ref_id).Cxx(:) = 0;
                else
                    Logger.getInstance.addError('Reference is missing, loosing all the coordinates');
                    for c = setdiff(1 : numel(coo_list), new_ref_id)
                        coo_list(c).rem(1:size(coo_list(c).xyz,1));
                    end
                end
            end
        end
        
        function computeModel(coo_list, enu_base)
            % Compute the spline model
            %
            % INPUT
            %   enu_base   [1x3] Array containing the dimensions of the spolien base for East, North and Up
            %
            % SYNTAX
            %   coo_list.computeModel(enu_anime);
            for coo = coo_list
                coo.toLocal;
                coo.splineModel(enu_base(1), enu_base(2), enu_base(3));
                coo.toECEF;
            end
        end
        
        function emptyModel(coo_list)
            % Delete the actually stored model
            %
            % SYNTAX
            %   coo_list.emptyModel();
            for coo = coo_list
                coo.xyz_model = [];
            end
        end
        
        function [m_data, sync_time, data_sync, n_obs] = getMeanReflectionHeight(coo_list, flag_height, out_thr)
            % Compute mean refllector heigth from a list of coordinnates
            % containing reflector heights
            %
            % this does not take in considration the heigth of the receiver
            % not the best solution, but reflectometry data are used to be
            % computed independently from the estimation of the coordinates
            % (to be improved in the future)
            %
            % INPUT
            %   coo_list    list of coordinates containing rh
            %   out_thr     outlier threshold around median
            %
            % SYNTAX
            %   [m_data, sync_time] = coo_list.getMeanReflectionHeight(thr)
            
            % threshold level
            if nargin < 3 || isempty(out_thr)
                out_thr = 0.2;
            end
            if nargin < 2 || isempty(flag_height)
                flag_height = true;
            end
            n_min_obs = 5;

            % collect time
            clear tmp;
            for i = 1 : numel(coo_list)
                tmp(i) = coo_list(i).getReflectometryTime;
            end
            
            % sync times
            [sync_time, id_sync] = tmp.getSyncedTime;
            % init syncd data
            data_sync = nan(size(id_sync));
            n_obs = zeros(size(id_sync));
            std_obs = zeros(size(id_sync));
            % sync data
            for i = 1:numel(coo_list)
                coo = coo_list(i);
                [~, ~, ~, h_masl] = coo_list(i).getGeodetic;

                time_coo = coo.getTime;
                time_data = coo.getReflectometryTime.getEpoch(noNaN(id_sync(:,i)));
                try
                    if flag_height
                        if time_coo.length > 1
                            data_sync(~isnan(id_sync(:,i)),i) = interp1(time_coo.getMatlabTime, h_masl, time_data.getMatlabTime, 'spline') - coo.reflectometry.value(noNaN(id_sync(:,i)));
                        else
                            data_sync(~isnan(id_sync(:,i)),i) = h_masl - coo.reflectometry.value(noNaN(id_sync(:,i)));
                        end
                    else
                        data_sync(~isnan(id_sync(:,i)),i) = coo.reflectometry.value(noNaN(id_sync(:,i)));
                    end

                    % if there is no reflectometry solution this might fail
                    n_obs(~isnan(id_sync(:,i)),i) = coo_list(i).reflectometry.n_obs(noNaN(id_sync(:,i)));
                    data_sync(n_obs(:,i) < n_min_obs, i) = nan; 
                    std_obs(~isnan(id_sync(:,i)),i) = coo_list(i).reflectometry.std(noNaN(id_sync(:,i)));
                catch
                    % it's ok no reflectometry data
                end
            end
            data_sync_bk = data_sync;
            id_sync(isnan(data_sync)) = NaN; % coordinates might exist but without reflectometry
            
            % select the column with more data to compute inter solution biases
            [~, id_best] = min(sum(isnan(id_sync)));
            
            % compute receiver biases (I don't have precise enough heigths to use them)
            biases = median(data_sync - data_sync(:, id_best), 1, 'omitnan');
            id_ok = ~isnan(biases);
            mean_bias = sum(biases(id_ok) .* sum(nan2zero(n_obs(:,id_ok)),1)) / sum(serialize(nan2zero(n_obs(:,id_ok))),1);
            biases = biases - mean_bias;
            
            % remove biases
            data_sync = bsxfun(@minus, data_sync_bk, biases);
            
            % outlier detection
            if sync_time.getRate < 43200
                w_obs = zero2nan(n_obs); % Estimate the weights on the base of the number of data and std dev
            else
                w_obs = (zero2nan(n_obs) ./ zero2nan(std_obs+1e-3).^2); % Estimate the weights on the base of the number of data and std dev
            end
            w_obs = nan2zero(w_obs ./ repmat(sum(nan2zero(w_obs),2), 1, size(w_obs,2))); % normalize weigths
            m_data = sum(nan2zero(data_sync) .* nan2zero(w_obs), 2); % compute weighted mean
            
            % outlier detection sensor
            sensor = abs(data_sync - m_data);
            thr = zeros(size(sensor,1),1); % compute a variable thr I always want at least one obs
            flag_ko = false(numel(thr),1);
            thr_lim = min(3*perc(sensor(:),0.8), out_thr);
            for i = 1 : numel(thr)
                if any(data_sync(i,:))
                    thr(i) = max(0.1, max(min(noNaN(sensor(i,:))), thr_lim));
                else
                    flag_ko(i) = true;
                end
            end
            sensor = sensor - repmat(thr, 1, size(sensor,2));
            id_ko = isnan(data_sync) | sensor > 0; % remove what I think are the outliers
            n_obs_new = n_obs;
            n_obs_new(id_ko) = 0;
            
            % Recompute mean without this suspected outliers
            biases = biases + median(data_sync - m_data, 1, 'omitnan');
            id_ok = ~isnan(biases);
            mean_bias = sum(biases(id_ok) .* sum(nan2zero(n_obs_new(:,id_ok)),1)) / sum(serialize(nan2zero(n_obs_new(:,id_ok))),1);
            biases = biases - mean_bias;
            
            % remove biases
            data_sync = bsxfun(@minus, data_sync_bk, biases);
            
            % weight of observations for the suspected outliers are at zero
            w_obs(id_ko) = 0;
            w_obs = nan2zero(w_obs ./ repmat(sum(nan2zero(w_obs),2), 1, size(w_obs,2))); % normalize weigths
            m_data = sum(nan2zero(data_sync) .* nan2zero(w_obs), 2);
            sensor = abs(data_sync - m_data);
            id_ko0 = id_ko;
            id_ko = isnan(data_sync) | sensor > out_thr | n_obs < n_min_obs; % remove what you think are the outliers
            
            % if all are considered outliers
            id_all = find(sum(id_ko,2) == size(id_ko, 2));
            id_ko(id_all,:) = id_ko0(id_all,:);
            if any(id_all)
                for i = 1:numel(id_all)
                    [~, id_ok] = max(n_obs(id_all(i),:));
                    id_ok = id_ok & n_obs(id_all(i),:) >= n_min_obs;
                    % Keep the best solution as the one computed with
                    % the greter number of observations
                    for r = find(id_ok)
                        id_ko(id_all(i) + (r-1) * size(id_ko,1)) = false;
                    end
                end
            end
            % filter out the outliers
            data_sync(id_ko) = nan;
            n_obs(id_ko) = 0;
            
            % mean weighted by the number of observations used to compute it
            if sync_time.getRate < 43200
                w_obs = zero2nan(n_obs.^2); % Estimate the weights on the base of the number of data and std dev
            else
                w_obs = (zero2nan(n_obs.^2) ./ zero2nan(std_obs+1e-3).^2);  % Estimate the weights on the base of the number of data and std dev
            end            
            w_obs = nan2zero(w_obs ./ repmat(sum(nan2zero(w_obs),2), 1, size(w_obs,2)));
            m_data = sum(nan2zero(data_sync) .* nan2zero(w_obs), 2);
            m_data(flag_ko | ~any(data_sync,2)) = nan;
            data_sync = data_sync_bk;
        end
        
        function [mask, fh] = getSkyMask(coo, varargin)
            % Get a predicted sky mask from DTM
            % Warning it does not handle yet areas close to 180 -180 longitudes
            %
            % SYNTAX
            %   [mask, fh] = coo.getSkyMask(...)
            %
            % INPUT
            %   Input are defined as a list of (couples option, value)
            %   e.g. coo.getSkyMask('rec', rec, 'flag:_fig', true, ...)
            %
            %   'rec'          use an external receiver
            %                      GNSS_Station
            %                      default:  empty
            %
            %   'use_dtm_h0'   use the height of the DTM as height of the
            %                  point, if folse use the height in rec/coo
            %                      flag:     true / false
            %                      default:  true
            %
            %   'h_antenna'    additional heigth of the antenna
            %                      numeric
            %                      default:  0
            %
            %   'h_vegetation'  additional heigth of the DTM
            %                      numeric
            %                      default:  0
            %
            %   'az_mask'      array of azimuth for the mask
            %                      numeric array [deg]
            %                      default:  (0:0.5:360)'
            %
            %   'az_res'       when az_mask is of one element use this number as resolution5295
            %                      numeric [deg]
            %                      default:  0.5
            %
            %   'min_el'       minimum elevation to display [deg]
            %                      numeric [deg]
            %                      default:  -7.5
            %
            %  - FIGURE ---------------------------------------------------
            %
            %   'show_fig'     display the results on a figure
            %                      flag:     true / false
            %                      default:  false
            %
            %   'name'         display a different name of the station
            %                      char
            %                      default:  name of coordinates
            %
            %   'marker'       display a different marker of the station on the map
            %                      char
            %                      default:  marker of coordinates
            %
            %  - ORBITS ---------------------------------------------------
            %    Note: orbits are computed using final GFZ products 
            %    (or the loaded core if a rec is passed)
            %
            %   'date_start'   first epoch to compute orbits
            %                      GPS_Time
            %                      default:  GPS_Time('2021-11-1 00:00:00')
            %
            %   'date_stop'    last epoch to compute orbits
            %                      GPS_Time
            %                      default:  GPS_Time('2021-11-1 23:59:59')
            % 
            %   'sys_list'     array of constellation to use
            %                      logical [7x1] || char
            %                      default:  'GRE'
            % 
            
            % Default args --------------------------------------------
            
            rec = [];
            flag_rec = false;
            flag_linear_interp = true;
            
            % Additional antenna heigth
            h_antenna = 1;
            h_vegetation = 0;
            flag_h_tdm = true;
            
            % flag for displaying the figure;
            flag_fig = false;
            flag_only_sun = false; % show opnly the sun in the terrain
            
            % Coordinate name            
            name = coo.getName;
            marker = coo.getNameV3;
            
            % Azimuth resolution
            az_mask = (0:0.5:360)';
                        
            % Satellite system to plot as logical array
            % GPS GLO GLA QZS BDS IRN SBS
            sys_list = [1 1 1 0 0 0 0];
            flag_change_sys = false;
            
            % DTM map margin
            m_x = 0.5;
            m_y = 0.5;
            
            % size of the scatter points
            point_size = 10;
            
            % Orbit date interval
            date_start = [];
            date_stop = [];
                        
            % Fun params ----------------------------------------------
            
            min_el = -7.5;
            
            % Parse args ----------------------------------------------
            
            args = varargin;
            if iscell(args) && numel(args) > 0
                if iscell(args{1})
                    args = args{1};
                end
            end
            if nargin > 1
                a = 1;
                while a < numel(args)
                    if ischar(args{a})
                        switch args{a}
                            case {'rec'}
                                a = a + 1;
                                rec = args{a};
                                coo = rec.getCoo;
                                flag_rec = true;
                            case {'use_dtm_h0'}
                                    a = a + 1;
                                    flag_h_tdm = logical(args{a});
                            case {'only_sun'}
                                    a = a + 1;
                                    flag_only_sun = logical(args{a});                                   
                            case {'h_antenna'}
                                    a = a + 1;
                                    h_antenna = args{a};
                            case {'h_vegetation'}
                                    a = a + 1;
                                    h_vegetation = args{a};
                            case {'min_el'}
                                    a = a + 1;
                                    min_el = args{a};
                            case {'flag_fig', 'show_fig'}
                                    a = a + 1;
                                    flag_fig = logical(args{a});
                            case {'name'}
                                    a = a + 1;
                                    name = char(args{a});
                            case {'marker'}
                                    a = a + 1;
                                    marker = char(args{a});
                            case {'az_mask', 'az_res'}
                                    a = a + 1;
                                    az_mask = args{a};
                                    if numel(az_mask) == 1
                                        az_mask = (0:az_mask:360)';
                                    end
                            case {'date_start'}
                                    a = a + 1;
                                    date_start = args{a};
                                    if isnumeric(date_start) || ischar(date_start)
                                        date_start = GPS_Time(date_start);
                                    end
                            case {'date_stop'}
                                    a = a + 1;
                                    date_stop = args{a};
                                    if isnumeric(date_stop) || ischar(date_stop)
                                        date_stop = GPS_Time(date_stop);
                                    end
                            case {'sys_list'}
                                a = a + 1;
                                if ischar(args{a})
                                    sys_list = false(1,7);
                                    for ss = args{a}
                                        switch ss
                                            case 'G'
                                                sys_list(1) = true;
                                            case 'R'
                                                sys_list(2) = true;
                                            case 'E'
                                                sys_list(3) = true;
                                            case 'Q'
                                                sys_list(4) = true;
                                            case 'C'
                                                sys_list(5) = true;
                                            case 'I'
                                                sys_list(6) = true;
                                            case 'S'
                                                sys_list(7) = true;
                                        end
                                    end
                                else
                                    sys_list = logical(args{a});
                                end
                                flag_change_sys = true;
                            otherwise
                                if ischar(args{a})
                                    log = Core.getLogger;
                                    log.addWarning(sprintf('Parameter "%s" unrecognized in Coordinates.showMap()', args{a}));
                                end
                                a = a + 1;
                        end
                    else
                        % Ignore parameter
                    end
                    a = a + 1;
                end
            end
                        
            if isempty(date_start)
                if not(isempty(rec))
                    try
                        if rec.getTime.isEmpty
                            time = rec.work.time.first;
                        else
                            time = rec.getTime.first;
                        end
                        date_start = time;
                    catch
                        date_start = GPS_Time('2022-04-20 00:00');
                    end
                else
                    date_start = GPS_Time('2022-04-20 00:00');
                end
            end
            if isempty(date_stop)
                if not(isempty(rec))
                    try
                        if rec.getTime.isEmpty
                            time = rec.work.time.last;
                        else
                            time = rec.getTime.last;
                        end
                        date_stop = time;
                    catch
                        date_stop = GPS_Time('2022-04-20 23:59:59');
                    end
                else
                    date_stop = GPS_Time('2022-04-20 23:59:59');
                end
            end
            
            log = Core.getLogger;
            coo = coo(not(coo.isEmpty));
            % Init output
            el_mask = zeros(size(az_mask));
            az_mask_out = zeros(size(az_mask));
            fh = [];
            
            if not(isempty(coo))
                if isempty(rec)
                    [lat, lon, h_ell] = coo.getMedianPos.getGeodetic();
                else
                    [lat, lon, h_ell] = rec.getPos.getMedianPos.getGeodetic();
                end
                lat = lat * 180/pi;
                lon = lon * 180/pi;
                
                % Get the DTM
                log.addMarkedMessage('Get DTM');
                % res = '1' / '3' / 'maximum' / 'high' / 'medium' / 'low'
                nwse = [min(90,lat + m_y) max(-180, lon - m_x) max(-90, lat - m_y) min(180, lon + m_x)];
                [dtm, lat_dtm, lon_dtm] = Core.getRefDTM(nwse, 'ortho', '3', false);
                dtm = fillmissing(dtm, 'linear'); % <= needed for griddedInterpolant
                % step of the DTM map
                step_lon = median(diff(lon_dtm));
                step_lat = median(diff(lat_dtm));
                %dtm = flipud(dtm);
                
                % Extract pers azimuth profile
                % [id_x, id_y] id of the center of the DTM map
                
                [~,idy] = min(abs(lat_dtm - lat));
                if numel(idy) > 0
                    idy = idy(1); % take only one point
                end
                
                [~,idx] = min(abs(lon_dtm - lon));
                if numel(idx) > 0
                    idx = idx(1); % take only one point
                end
                
                % Height of the point on the DTM (nearest neighbour)
                h0_dtm = nan2zero(double(dtm(idy, idx)));
                
                % Height of the point on the DTM (interp using close points)
                sub_idy = idy + [-5:5]; sub_idy(sub_idy < 1 | sub_idy > size(dtm,1)) = [];
                sub_idx = idx + [-5:5]; sub_idx(sub_idx < 1 | sub_idx > size(dtm,2)) = [];
                x_grid = lon_dtm(sub_idx);
                y_grid = lat_dtm(sub_idy);
                [xmg, ymg] = meshgrid(x_grid, y_grid);
                finterp = scatteredInterpolant(xmg(:), ymg(:), double(serialize(dtm(sub_idy, sub_idx))), 'natural');
                h0_dtm = finterp(lon, lat) + h_antenna;
                
                % If a receiver is passed let's use it's coordinates
                if ~flag_h_tdm
                    h0 = h_ell + h_antenna;
                    log.addMarkedMessage(sprintf('Considering ellipsoidal point heigth to %.1fm, the DTM say it should be %.1fm', h_ell, h0));
                else
                    h0 = h0_dtm + h_antenna;
                    log.addMarkedMessage(sprintf('Considering ellipsoidal point heigth to %.1fm from the DTM, the coordinates say it should be %.1fm', h0, h_ell));
                end
                
                % Create the coordinates of the origin point
                coo_origin = Coordinates.fromGeodetic(lat/180*pi, lon/180*pi, h0);
                coo_origin.setName(marker);
                
                if flag_fig
                    % Prepare the figure;
                    fh = figure('Visible', 'off');
                    Core_UI.resizeFig(fh, [1600,900]);
                    
                    ax_map = subplot(4,3,[1 4]); axis;
                    ax_sky = subplot(4,3,[2 5]); axis;
                    ax_sky_tracks = subplot(4,3,[3 6]); axis;
                    ax_terrain = subplot(4,3,[7:12]); axis;
                    
                    % Prepare terrain data
                    az_terrain = [];
                    el_terrain = [];
                    d_terrain = [];
                    size_dot = [];
                end
                
                % Prepare elevation mask
                a = 0;
                % Compute true distance between the point and the coordinates in the direction of the azimuth
                log.addMarkedMessage(sprintf('Computing the mask @ %.2f degree of resolution', median(diff(az_mask))));
                
                % Create DTM gridded interpolant
                if flag_linear_interp
                    [Y,X] = ndgrid(lat_dtm,lon_dtm);
                    fdtm = griddedInterpolant(X',Y',dtm');
                end

                for az = az_mask'
                    a = a + 1;

                    if ~flag_linear_interp
                        % Intercept the knot of the DTM intersecting a line starting from the central node
                        % For each lat get a lon intersecting the line of sight
                        if (az >= 0 && az < 45) || (az >= 135 && az < 225) || (az >= 315 && az <= 360)
                            if (az < 90) || (az > 270)
                                idy_l = idy + 1 : numel(lat_dtm); % for all the lat above
                            else
                                idy_l = idy - (1 : (idy -1)); % for all the lat above
                            end
                            d_lon = (lat_dtm(idy_l) - lat) .* tand(az);
                            idx_l = idx + round(d_lon / step_lon);
                            if az < 180
                                id_ko = idx_l > numel(lon_dtm);
                            else
                                id_ko = idx_l <= 0;
                            end
                            % For each lon get a lat intersecting the line of sight
                        elseif (az >= 45 && az < 135) ||  (az >= 225 && az < 315)
                            if (az >= 45 && az < 135)
                                idx_l = idx + 1 : numel(lon_dtm); % for all the lat above
                            else
                                idx_l = idx - (1 : (idx -1)); % for all the lat above
                            end
                            d_lat = (lon_dtm(idx_l) - lon) .* cotd(az);
                            idy_l = idy + round(d_lat / step_lat);
                            if az < 90 || az > 270
                                id_ko = idy_l > numel(lat_dtm);
                            else
                                id_ko = idy_l <= 0;
                            end
                        end

                        % Shorten the line of sight if it goes out of the map
                        idy_l(id_ko) = [];
                        idx_l(id_ko) = [];

                        % coordinates on the line of sight

                        h_list = nan2zero(double(dtm(idy_l(:) + (idx_l(:)-1) * size(dtm,1))));
                        r1 = h0 + GPS_SS.R_EARTH;
                        r2 = nan2zero(double(dtm(idy_l(:) + (idx_l(:)-1) * size(dtm,1)))) + GPS_SS.R_EARTH;

                    else % linear interp

                        % Line of sight positions
                        step = mean([diff(lat_dtm); diff(lon_dtm)])/2;

                        max_len = sqrt(max(abs(minMax(lat_dtm-lat))).^2 + max(abs(minMax(lon_dtm-lon))).^2) / step;
                        lat_line = ((1:max_len) * cosd(az)) * step + lat;
                        lon_line = ((1:max_len) * sind(az)) * step + lon;
                        id_ok = lat_line(:) >= min(lat_dtm) & ...
                            lat_line(:) <= max(lat_dtm) & ...
                            lon_line(:) >= min(lon_dtm) & ...
                            lon_line(:) <= max(lon_dtm);
                        lat_line = lat_line(id_ok)';
                        lon_line = lon_line(id_ok)';
                        line_val = fdtm(lon_line, lat_line);

                        % coordinates on the line of sight                   

                        % reduce to the nuber of point maybe visible
                        h_list = nan(size(line_val));
                        h_list(1) = line_val(1);
                        for i = 2:numel(line_val)
                            h_list(i) = max(h_list(i-1), line_val(i));
                        end
                        id_ko = (line_val > h0) & [false; diff(h_list) == 0];
                        h_list(line_val <= h0) = line_val(line_val <= h0);
                        
                        h_list(id_ko) = nan;
                        lat_line(id_ko) = nan;
                        lon_line(id_ko) = nan;                        
                        r1 = h0 + GPS_SS.R_EARTH;
                        r2 = line_val + GPS_SS.R_EARTH;
                    end

                    % use geometrical distance
                    if flag_linear_interp
                        [el, d, true_az] = Coordinates.getSigthElevation(lat, lon, r1, lat_line, lon_line, r2 + h_vegetation);
                    else
                        [el, d, true_az] = Coordinates.getSigthElevation(lat, lon, r1, lat_dtm(idy_l), lon_dtm(idx_l), r2 + h_vegetation);
                    end
                    el = max(min_el, el);
                    
                    % Compute elevation mask
                    if ~isempty(el)
                        [el_mask(a), id_max] = max(noNaN(el));
                        az_mask_out(a) = true_az(id_max);
                        %az_mask_out(a) = median(true_az, 'omitnan');
                        
                        if flag_fig
                            % Draw the mountains
                            % Keep for each elevation only the points that can be seen
                            setAxis(fh,4);
                            el_max = min_el;
                            id_ok = ~isnan(el);
                            for i = find(id_ok(:)')
                                id_ok(i) = el_max < el(i);
                                if id_ok(i)
                                    el_max = el(i);
                                end
                            end
                            %az_terrain = [az_terrain; repmat(az,sum(id_ok),1)];
                            %az_terrain = [az_terrain; repmat(median(true_az),sum(id_ok),1)];
                            az_terrain = [az_terrain; true_az(id_ok)];
                            el_terrain = [el_terrain; el(id_ok)];
                            d_terrain = [d_terrain; d(id_ok)];
                            size_dot = [size_dot; diff([find(id_ok); numel(id_ok)])];
                        end
                    end
                end

                % Reorder and reallign mask
                az_mask_out(end) = []; % this is 360 == 0
                el_mask(end) = [];     % this is 360 == 0                
                [az_sort, ids] = sort(az_mask_out);
                el_sort = el_mask(ids);
                el_mask = interp1q([az_sort(end-49:end) - 360; az_sort; az_sort(1:50) + 360], ...
                         [el_sort(end-49:end); el_sort; el_sort(1:50)], ...
                         az_mask);  
                
                if flag_fig
                    % Set the mountains color
                    setAxis(fh,4);
                    cmap = flipud(Cmap.get('PuBuGn', 1400)); cmap = cmap(1:1024,:); colormap(cmap); colorbar;

                    % Plot sun path ----------------------------------------------------------------------
                    if sum(sys_list) == 0 || flag_only_sun
                        ax_cartesian = ax_terrain;
                    else
                        ax_cartesian = [];
                    end

                    % Sun in different time of the year --------------------------------------------------

                    % Equinoxium
                    time_sun = GPS_Time(datenum('2024-09-23 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, ax_sky, true, false);

                    % Solstitium
                    time_sun = GPS_Time(datenum('2021-12-21 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, ax_sky, true, lat <= 0);

                    time_sun = GPS_Time(datenum('2022-06-21 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, ax_sky, true, lat > 0);

                    % Others
                    time_sun = GPS_Time(datenum('2022-07-21 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, [], true, false);

                    time_sun = GPS_Time(datenum('2022-08-21 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, [], true, false);

                    time_sun = GPS_Time(datenum('2022-10-21 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, [], true, false);

                    time_sun = GPS_Time(datenum('2022-11-21 00:00')+[0:1/288:1]');
                    coo.sunPlot(time_sun, min_el, ax_cartesian, [], true, false);
                    
                    subplot(ax_terrain);

                    % Plot terrain -----------------------------------------------------------------------
                    
                    % plot the terrain
                    p = patch([az_mask; fliplr(minMax(az_mask))'], [el_mask; min_el; min_el], cmap(1,:));
                    p.FaceAlpha = 0.6;

                    % fh.Visible = 'on'; subplot(ax_terrain); plot(az_mask, el_mask, '-k', 'LineWidth', 2);
                    
                    % set the style of the figure
                    ax_terrain.GridAlpha = 0.11;
                    ax_terrain.Box = 'on';

                    xlabel('Azimuth [deg]');
                    ylabel('Elevation [deg]');
                    title(sprintf('%s  -  %.6f, %.6f @ %.1f m (DTM %.1fm)', name, lat, lon, h0, h0_dtm));
                    cb = colorbar; title(cb, sprintf('Distance\n[Km]'));
                    ax_terrain.XTick = 0:5:360;
                    ax_terrain.XTickLabel = {'0 (N)', '', '', '', '', '', '', '', '', ...
                        '45 (NE)', '', '', '', '', '', '', '', '', ...
                        '90 (E)', '', '', '', '', '', '', '', '', ...
                        '135 (SE)', '', '', '', '', '', '', '', '', ...
                        '180 (S)', '', '', '', '', '', '', '', '', ...
                        '225 (SW)', '', '', '', '', '', '', '', '', ...
                        '270 (W)', '', '', '', '', '', '', '', '', ...
                        '315 (NW)', '', '', '', '', '', '', '', '', ...
                        '360 (N)'};
                    ax_terrain.YTick = 0:5:90;
                    ax_terrain.YTickLabel = {'0', '', '10', '', '20', '', '30', '', '40', '', '50', '', '60', '', '70', '', '80', '', '90'};
                    ax_terrain.LineWidth = 1.5;
                    ax_terrain.Box = 'on';
                    ax_terrain.XGrid = 'on';
                    ax_terrain.YGrid = 'on';
                    ax_terrain.XAxis.TickLength = [0.002 0.0001];
                    ax_terrain.YAxis.TickLength = [0.002 0.0001];
                    ylim([min_el max(90, max(el_mask))]);
                    xlim([0 360]);
                    
                    % plot terrain
                    point_size = min(6,size_dot);
                    point_size = round((max(0,(point_size-1))/max(point_size-1)) * 10)+0.85;
                    hold on; scatter(az_terrain, el_terrain, point_size, min(max(10,perc(d_terrain,0.97)),d_terrain)/1e3, 'filled'); hold on;
                    
                    % Plot horizon -----------------------------------------------------------------------
                    
                    plot([0 360], [0 0], '-.', 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;      
                    drawnow
                end
                
                % set elevation cut-off to zero, this is done since min_el could be < 0
                el_mask(a) = max(0, el_mask(a));
                if flag_fig
                    % display the predicted mask
                    setAxis(fh,2);
                    if ~isempty(rec)
                        [snr_map, snr_map_fill, snr_mask, n_data_map, out_map] = GReD_Utility.getMeanMapSNR(rec, 0.5, 32, 22, unique(rec.work.system));
                        if any(snr_mask(:))
                            [~, id_min] = max(snr_mask); id_min(id_min == 1) = 180; id_min = 90-(id_min/2);
                            az_grid = linspace(0,360,720);
                            p = polarPatch(az_grid / 180 * pi, (90 - id_min) / 180 * pi, [0.5 0.5 0.5], false); p.FaceAlpha = 0.4; 
                            hold on;
                        end
                    end
                    if isempty(sys_list)
                        sys_list = [1 1 1 0 0 0 0];
                    end
                    if sys_list(2)
                        ss = GLONASS_SS;
                    elseif sys_list(1)
                        ss = GPS_SS;
                    elseif sys_list(5)
                        ss = BeiDou_SS;
                    elseif sys_list(3)
                        ss = Galileo_SS;
                    elseif sys_list(4)
                        ss = QZSS_SS;
                    elseif sys_list(6)
                        ss = IRNSS_SS;
                    else
                        ss = GPS_SS;
                    end
                        
                    % Polar mask due to satellite systems orbit inclinations
                    [mn, ms, msk] = ss.getPolarMask(coo(1).lat, coo(1).lon, 0, 0.5);
                    p = polarPatch(az_mask / 180 * pi, (90-el_mask) / 180 * pi, Core_UI.LBLUE, false);
                    p.FaceAlpha = 0.5;
                    p.EdgeAlpha = 0.5;
                    color = [0.6 0.2 0.2];
                    color_border = [0.8 0.1 0.1];
                    if ~isempty(mn)
                        p = polarPatch(mn(:,1), pi/2-mn(:,2), color, isempty(ms)); hold on;
                        p.FaceAlpha = 0.5;
                        p.EdgeColor = color_border;
                        p.EdgeAlpha = 0.5;
                    end
                    if ~isempty(ms)
                        p = polarPatch(ms(:,1), pi/2-ms(:,2), color, true); hold on;
                        p.FaceAlpha = 0.5;
                        p.EdgeColor = color_border;
                        p.EdgeAlpha = 0.5;
                    end
                    title(sprintf('Predicted GNSS satellite visibility\\fontsize{5} \n'));
                end
                
                % ORBIT DETERMINANTION
                if isempty(sys_list)
                    sys_list = [1 1 1 0 0 0 0];
                end
                        
                if flag_fig
                    if any(sys_list)
                        log.addMarkedMessage('Orbit determination');

                        core = Core.getCurrentCore;
                        cc = Core.getConstellationCollector();

                        if isempty(rec)
                            cc.setActive(sys_list);

                            % Create a fake receiver
                            rec = GNSS_Station();
                            core.rec = rec;
                            % set its coordinates
                            rec.work.xyz = coo_origin.getXYZ;
                            % Set its name
                            rec.setMarkerName(name);
                            % Fake observation time, one per minute
                            rec.work.time = GPS_Time.fromRefTime(date_start.getMatlabTime, 0:60:(date_stop - date_start));

                            n_epoch = rec.work.time.length;
                            n_sat = numel(cc.index);

                            % Set fake observations (SNR) => needed later to compute the orbits
                            rec.work.go_id = cc.index;
                            rec.work.obs_code = repmat('S  ', n_sat, 1);
                            rec.work.obs = ones(n_sat, n_epoch);

                            % Prepare the computational enviroment (BREVA core)
                            state = core.getState;
                            % Set session limits
                            state.setSessionStart(date_start);
                            state.setSessionStop(date_stop);
                            state.setSessionDuration(date_stop-date_start);
                            state.setBuffer(0,0);
                            state.setRinexSession(false);
                            % Set the orbit type e.g. GFZ final
                            %state.setPreferredOrbit(5, 'code_predicted5')
                            %state.setPreferredOrbit(1, 'gfz'); % broadcast
                            %state.setPreferredOrbit(1, 'code_mgex_aiub_rnx3'); % broadcast
                            state.setPreferredOrbit(1, 'mgex_broadcast'); % broadcast
                            core.sky.clearOrbit();
                            state.updateNavFileName();
                            % Download the orbits if needed
                            fw = File_Wizard;
                            fw.conjureNavFiles(date_start, date_stop);
                            % Prepare the Core Sky object (load the orbits)
                            state.setCurSession(1);
                            core.sky.initSession(date_start, date_stop, cc, true);
                        else
                            if ~flag_change_sys
                                sys_list = cc.getActive;
                            end
                        end

                        % Update AZ and Elevation
                        if not(any(rec.work.sat.az(:)))
                            rec.work.updateAzimuthElevation();
                        end

                        [rec_az, rec_el] = rec.work.getAzEl;
                        rec_az = mod(rec_az, 360);

                        % filter az, el => apply mask
                        id_ko = isnan(zero2nan(rec_az));
                        if not(flag_rec)
                            az_step = median(diff(az_mask));
                            for a = 1 : numel(az_mask)
                                id_ko = id_ko | round((rec_az / az_step)) == round(az_mask(a) / az_step) & rec_el < el_mask(a);
                            end
                        end
                        id_ko(:,ismember(cc.system, cc.SYS_C(~sys_list))) = false; % remove constellations not enabled / requested
                        rec_az(id_ko) = nan;
                        rec_el(id_ko) = nan;

                        % compute constellation color
                        ssys = zeros(size(rec_az));
                        ssys(:, cc.system == 'G') = 1;
                        ssys(:, cc.system == 'R') = 2;
                        ssys(:, cc.system == 'E') = 3;
                        ssys(:, cc.system == 'J') = 5;
                        ssys(:, cc.system == 'C') = 4;
                        ssys(:, cc.system == 'I') = 6;
                        ssys(:, cc.system == 'S') = 7;

                        col = Core_UI.getColor(1:7,7);

                        % Plot Satellite orbits
                        if ~flag_only_sun
                            ax = setAxis(fh,4);
                            scatter(rec_az(~id_ko), rec_el(~id_ko), 15, col(ssys(~id_ko),:), 'filled');
                        end

                        % Polar plot orbits
                        ax = setAxis(fh,3);
                        polarScatter(rec_az(~id_ko)./180*pi, (90-rec_el(~id_ko))./180*pi, 15, col(ssys(~id_ko),:), 'filled');
                        title('Daily tracks');
                    else % No orbits to plot
                        ax = setAxis(fh,3);
                        polarScatter(0, 0, 0.001, [0 0 0]);
                        title('No constellations selected');
                    end
                end
                
                if flag_fig
                    log.addMarkedMessage('Add the map of the point');
                    % Add the map
                    setAxis(fh,1);
                    ax = setAxis(fh,1); coo_origin.showMap('new_fig', false, 'bg_type', 'map', 'add_margin', [-0.029 -0.046], 'contour_lines', false);
                    % Move patch down:
                    patches = false(numel(ax.Children),1); 
                    for c = 1:numel(ax.Children)
                        patches(c) = isa(ax.Children(c), 'matlab.graphics.primitive.Patch'); 
                    end
                    if any(patches)
                        ax.Children(find(patches,1,'first')).Vertices(:,2) = ax.Children(find(patches,1,'first')).Vertices(:,2) - 0.000005;
                        for i = 4:find(patches,1,'first')-1
                            ax.Children(i).Position(2) = ax.Children(i).Position(2) - 0.000005;
                        end
                    end
                    
                    % Set figure name
                    fh.Name = sprintf('%03d: Mask %s', fh.Number, strrep(name,' ', '_')); fh.NumberTitle = 'off';
                    fig_name = sprintf('PredictedMask_%s', name);
                    fh.UserData = struct('fig_name', fig_name);
                    
                    % Beautify figure
                    Core_UI.insertLogo(fh);
                    Core_UI.beautifyFig(fh);
                    Core_UI.addBeautifyMenu(fh);
                    Core_UI.addExportMenu(fh);
                    fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                end
                log.addMarkedMessage('The elevation mask is ready');
            end
            
            mask = [az_mask(:), el_mask(:)];
        end   
    end

    methods (Access = 'private')
        function hl = sunPlot(coo, time_sun, min_el, cartesian_ax, polar_ax, flag_points, flag_labels)
            % Plot the sun path
            %
            % SYNTAX
            %   hl = sunPlot(time_sun)
            sky = Core.getCoreSky;
            [sun_el, sun_az] = sky.sunElevation(repmat(coo.xyz,time_sun.length(),1), time_sun);
            id_ko = (sun_el < min_el);
            sun_el(id_ko) = nan; sun_az(id_ko) = nan;
            if ~isempty(polar_ax)
                subplot(polar_ax);
                hl = polarPlot(mod(sun_az, 360)./180*pi, (90-sun_el)./180*pi, '-', 'LineWidth', 1.4, 'Color', Core_UI.YELLOW); hold on;                
            end
            if ~isempty(cartesian_ax)
                subplot(cartesian_ax);
                hl = plot(mod(sun_az, 360), sun_el, '-', 'LineWidth', 2, 'Color', Core_UI.YELLOW); hold on;
            end
            if flag_points
                % Add a point for each hours
                id_hours = find(mod(round(time_sun.getMatlabTime*24e3)/1e3,1) == 0);
                az = mod(sun_az(id_hours), 360);
                el = sun_el(id_hours);
                if ~isempty(polar_ax)
                    subplot(polar_ax);
                    polarPlot(mod(az, 360)./180*pi, (90-el)./180*pi, '.', 'MarkerSize', 20, 'Color', Core_UI.YELLOW);
                    for j = 1:numel(el)
                        if el(j) > 0
                            decl_n = ((90-el(j))/180*pi)/(pi/2);
                            x = sind(az(j)) .* decl_n;
                            y = cosd(az(j)) .* decl_n;
                            text(x,y + 0.1, time_sun.getEpoch(id_hours(j)).toString('HH'), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [224/256 186/256 82/256], 'FontWeight', 'bold');
                        end
                    end
                end
                if ~isempty(cartesian_ax)
                    subplot(cartesian_ax);
                    plot(az, el, '.', 'MarkerSize', 25, 'Color', Core_UI.YELLOW);
                    if flag_labels
                        for j = 1:numel(el)
                            if el(j) > min_el
                                text(az(j),el(j) + 3.5, time_sun.getEpoch(id_hours(j)).toString('HH'), 'HorizontalAlignment', 'center', 'FontSize', 9, 'Color', [224/256 186/256 82/256], 'FontWeight', 'bold');
                            end
                        end
                    end
                end
            end
        end
    end
    
    % =========================================================================
    %%    COORDINATES OPERATIONS
    % =========================================================================
    
    methods (Access = 'public', Static)

        function distances = vincentyDistance(lat_ref, lon_ref, lat_list, lon_list)
            % Distance on the ellipsoid in km
            %
            % INPUT
            %   lat_ref     latitude of the reference point
            %   lon_ref     longitude of the reference point
            %   lat_list    latitude of a list of point to compute distances with
            %   lon_ref     longitude of a list of point to compute distances with
            %
            % OUTPUT
            %   distances   distances in kilometers
            %
            % SYNTAX
            %   distances = Coordinates.vincentyDistance(lat_ref, lon_ref, lat_list, lon_list)

            % WGS84 ellipsoid parameters
            a = 6378137; % semi-major axis in meters
            f = 1/298.257223563; % flattening
            b = (1-f)*a; % semi-minor axis

            % Convert degrees to radians
            lat_ref = deg2rad(lat_ref);
            lon_ref = deg2rad(lon_ref);
            lat_list = deg2rad(lat_list);
            lon_list = deg2rad(lon_list);

            U1 = atan((1-f) * tan(lat_ref));
            U2s = atan((1-f) * tan(lat_list));
            Ls = lon_list - lon_ref;
            Lambdas = Ls;
            sinU1 = sin(U1);
            cosU1 = cos(U1);
            sinU2s = sin(U2s);
            cosU2s = cos(U2s);

            iterLimit = 1000;

            for iteration = 1:iterLimit
                sinLambdas = sin(Lambdas);
                cosLambdas = cos(Lambdas);
                sinSigmas = sqrt((cosU2s.*sinLambdas).^2 + (cosU1*sinU2s-sinU1*cosU2s.*cosLambdas).^2);
                cosSigmas = sinU1*sinU2s + cosU1*cosU2s.*cosLambdas;
                sigmas = atan2(sinSigmas, cosSigmas);
                sinAlphas = cosU1 .* cosU2s .* sinLambdas ./ sinSigmas;
                cos2Alphas = 1 - sinAlphas.^2;
                cos2SigmasM = cosSigmas - 2*sinU1*sinU2s./cos2Alphas;
                C = f/16*cos2Alphas.*(4+f*(4-3*cos2Alphas));
                Lambdas_prev = Lambdas;
                Lambdas = Ls + (1-C) .* f .* sinAlphas .* (sigmas + C.*sinSigmas.*(cos2SigmasM+C.*cosSigmas.*(-1+2.*cos2SigmasM.^2)));
                if all(abs(Lambdas - Lambdas_prev) < 1e-12)
                    break;
                end
            end

            u2s = cos2Alphas * (a^2 - b^2) / (b^2);
            A = 1 + u2s/16384.*(4096+u2s.*(-768+u2s.*(320-175*u2s)));
            B = u2s/1024 .* (256+u2s.*(-128+u2s.*(74-47*u2s)));
            deltaSigmas = B.*sinSigmas.*(cos2SigmasM+B/4.*(cosSigmas.*(-1+2.*cos2SigmasM.^2)- B/6.*cos2SigmasM.*(-3+4.*sinSigmas.^2).*(-3+4.*cos2SigmasM.^2)));

            % Distance in meters
            s = b.*A.*(sigmas-deltaSigmas);

            % Convert to kilometers and return
            distances = s / 1000;
        end

    end
    % =========================================================================
    %%    COORDINATES TRANSFORMATIONS
    % =========================================================================
    
    methods (Access = 'public', Static)
        function pos_diff = rotateDisplacements(pos_diff, rot_angle_cw)
            % Input:
            %       pos_diff - n x 3 matrix containing the coordinates in ENU format
            %       rot_angle_cw - clockwise rotation angle in radians
            %  SINTAX:
            %       pos_diff = rotateDisplacements(pos_diff, rot_angle_cw)

            % Extract the East and North coordinates from pos_diff
            if ~isempty(pos_diff)
                east = pos_diff(:, 1);
                north = pos_diff(:, 2);

                % Compute the rotated coordinates
                cos_angle = cosd(rot_angle_cw);
                sin_angle = sind(rot_angle_cw);
                east_rotated = east * cos_angle + north * sin_angle;
                north_rotated = -east * sin_angle + north * cos_angle;

                % Update the rotated coordinates in pos_diff
                pos_diff(:, 1) = east_rotated;
                pos_diff(:, 2) = north_rotated;
            end
        end

        function N = getOrthometricCorrFromLatLon(phi, lam, geoid, method)
            % SYNTAX:
            %   N = getOrthometricCorr(phi, lam, geoid);
            %
            % EXAMPLE:
            %   core = Core.getInstance;
            %   core.initGeoid();
            %   Coordinates.getOrthometricCorrFromLatLon(45.69 ./ 180*pi, 9.03 ./ 180*pi)
            %   % answer should be 46.1767008
            %
            % INPUT:
            %   phi     = geodetic latitude                [deg (rad only for legacy method)]
            %   lam     = geodetic longitude               [deg (rad only for legacy method)]
            %   geoid   = regular map in geocentric coordinates <default = EGM08 0.5x0.5 deg>
            %   method  = interpolation approach:
            %              - legacy
            %              - grid
            %              - grid_cubic  ( phi, lam are array of grid coordinates )
            %              - grid_akima  ( phi, lam are array of grid coordinates )
            %              - linear
            %              - natural
            %
            % OUTPUT:
            %   N       = geoid ondulation [m]
            %
            % DESCRIPTION:
            %   Get the geoid ondulation (orthometric correction)
            
            if (nargin < 3) || isempty(geoid)
                geoid = Core.getRefGeoid();
            end
            
            if (geoid.grid == 0)
                core = Core.getInstance(false);
                core.initGeoid();
                geoid = core.getRefGeoid();
            end
            
            if (nargin < 4) || isempty(method)
                method = 'legacy';
            end
            
            if  (geoid.ncols == 0 || geoid.nrows == 0)
                Core.initGeoid();
                geoid = Core.getRefGeoid();
            end
            N = zeros(numel(lam), 1);
            
            switch method
                case 'legacy'
                    for i = 1 : numel(lam)
                        N(i) = grid_bilin_interp(lam(i) / pi * 180, phi(i) / pi * 180, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
                    end
                case 'grid'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
                    
                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    N = interp2(xmg, ymg, geoid.grid, lam, phi, 'linear');
                case 'grid_cubic'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
                    
                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    N = interp2(xmg, ymg, geoid.grid, lam, phi, 'cubic');
                case 'grid_akima'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
                    
                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    N = interp2(xmg, ymg, geoid.grid, lam, phi, 'makima');
                case 'linear'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    finterp = scatteredInterpolant(xmg(:), ymg(:), geoid.grid(:), 'linear');
                    
                    N = finterp(lam, phi);
                case 'natural'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
                    
                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    finterp = scatteredInterpolant(xmg(:), ymg(:), geoid.grid(:), 'natural');
                    N = finterp(lam, phi);
            end
        end
        
        function [lat, lon, h, lat_geoc] = cart2geod(xyz, y, z)
            % Get Geodetic coordinates from xyz cartesian coordinates
            %
            % SYNTAX:
            %   [lat, lon, h, lat_geoc] = Coordinates.cart2geod(xyz)
            %   [lat, lon, h, lat_geoc] = Coordinates.cart2geod(x, y, z)
            if isempty(xyz)
                lat = [];
                lon = [];
                h = [];
                lat_geoc = [];
            else
                if nargin == 1
                    z = xyz(:, 3);
                    y = xyz(:, 2);
                    xyz = xyz(:, 1);
                end
                
                a = Coordinates.ELL_A;
                e = Coordinates.ELL_E;
                
                % radius computation
                r = sqrt(xyz.^2 + y.^2 + z.^2);
                
                % longitude
                lon = atan2(y,xyz);
                
                % geocentric latitude
                lat_geoc = atan(z./sqrt(xyz.^2 + y.^2));
                
                % Coordinates transformation
                psi = atan(tan(lat_geoc)/sqrt(1-e^2));
                
                lat = atan((r.*sin(lat_geoc) + e^2*a/sqrt(1-e^2) * (sin(psi)).^3) ./ ...
                    (r.*cos(lat_geoc) - e^2*a * (cos(psi)).^3));
                
                if nargout > 2
                    N = a ./ sqrt(1 - e^2 * sin(lat).^2);
                    
                    % height
                    h = r .* cos(lat_geoc)./cos(lat) - N;
                end
            end
        end
        
        function [lat_geoc, lon, r] = cart2geoc(xyz, y, z)
            % Get Geocentric sphericals coordinates from xyz cartesian coordinates
            %
            % SYNTAX:
            %   [lat_geoc, lon, r] = Coordinates.cart2geoc(xyz)
            %   [lat_geoc, lon, r] = Coordinates.cart2geoc(x, y, z)
            
            if nargin == 1
                z = xyz(:, 3);
                y = xyz(:, 2);
                xyz = xyz(:, 1);
            end
            
            a = Coordinates.ELL_A;
            e = Coordinates.ELL_E;
            
            % radius computation
            r = sqrt(xyz.^2 + y.^2 + z.^2);
            
            % longitude
            lon = atan2(y,xyz);
            
            % geocentric latitude
            lat_geoc = atan(z./sqrt(xyz.^2 + y.^2));
        end
        
        function [east, north, utm_zone] = geod2utm(lat, lon)
            % Conversion from geodetic coordinates to planimetric coordinates (UTM WGS84).
            %
            % SYNTAX:
            %   [east, north, utm_zone] = geod2utm(lat, lon);
            %
            % INPUT:
            %   lat = latitude [rad]
            %   lon = longitude [rad]
            %
            % OUTPUT:
            %   east     = east Coordinates [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone
            %
            
            %number of input points
            n = size(lat);
            n = n(1,1);
            
            if (n > 0)
                % pre-allocation
                M = zeros(n,1);
                north_south = zeros(n,1);
                utm_zone(n,:) = '60 X';
                north = zeros(n,1);
                east = zeros(n,1);
                
                % conversion algorithm
                % UTM parameters
                EsMc = 500000; % false East
                
                % WGS84 ellipsoid parameters
                % semi-major equatorial axis [m]
                ell_a = Coordinates.ELL_A;
                
                % squared eccentricity (a^2-b^2)/a^2
                e2 = Coordinates.E2;
                
                % contraction factor
                contr = 0.9996;
                
                ell_a = ell_a * contr;
                ecc4 = e2 * e2;
                ecc6 = ecc4 * e2;
                ecc8 = ecc6 * e2;
                
                k0 = ell_a * (e2 / 4 + ecc4 * 3 / 64 + ecc6 * 5 / 256 + ecc8 * 175 / 16384);
                k = ell_a - k0;
                k1 = ell_a * (ecc4 * 13 / 96 + ecc6 * 59 / 384 + ecc8 * 1307 / 8192);
                k2 = ell_a * (ecc6 * 61 / 484 + ecc8 * 609 / 2048);
                k3 = ell_a * (ecc8 * 49561 / 322560);
                c1 = (e2 * 5 - ecc4) / 6;
                c2 = (ecc4 * 104 - ecc6 * 45) / 120;
                c3 = ecc6 * 1237 / 1260;
                
                % Sines, cosines and latitude powers
                Elix = [lat lon];
                latsessadec = lat ./ pi .* 180;
                lonsessadec = lon ./ pi .* 180;
                
                fiSin(:,1) = sin(Elix(:,1));
                fiCos(:,1) = cos(Elix(:,1));
                fiSin2(:,1) = fiSin .* fiSin;
                fiSin4(:,1) = fiSin2 .* fiSin2;
                fiSin6(:,1) = fiSin4 .* fiSin2;
                
                % UTM zone finding
                % https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system#/media/File:Universal_Transverse_Mercator_zones.svg
                M(:, 1) = fix((180 + lonsessadec(:, 1)) / 6) + 1;
                north_south(:,1) = latsessadec(:, 1) >= 0;
                UTM_LETTERS = 'CDEFGHJKLMNPQRSTUVWX';
                letter_id = uint8(floor((latsessadec + 80) / 8) + 1);
                pole_id = latsessadec > 84;
                letter_id(pole_id) = letter_id(pole_id) + uint8(1) + uint8(lonsessadec(pole_id) > 0);
                pole_id = latsessadec < -80;
                letter_id(pole_id) = letter_id(pole_id) - uint8(1) - uint8(lonsessadec(pole_id) <= 0);
                utm_zone(:,1:2) = num2str(nan2zero(M(:,1)), '%02d');
                utm_zone(letter_id > 0,4) = UTM_LETTERS(min(letter_id(letter_id > 0),20));
                utm_zone(letter_id == 0, :) = repmat('00 0', sum(letter_id == 0), 1);
                
                east(letter_id == 0) = nan;
                north(letter_id == 0) = nan;
                for i = find(letter_id > 0)'
                    % Longitude of the central meridian
                    LonMeridianoCentrale = -177 + 6 * (M(i, 1) - 1);
                    
                    % Distance of the point from the central meridian
                    % la_sd --> distance in decimal degrees
                    % la    --> distance in radians
                    la_sd = lonsessadec(i, 1) - LonMeridianoCentrale;
                    la = la_sd / 180 * pi;
                    
                    if la == 0
                        laSin = 0;
                        laCos = 1;
                        laCot = 0;
                    else
                        laSin = sin(la);
                        laCos = cos(la);
                        laCot = laCos / laSin ;
                    end
                    
                    % longitude with respect to central meridian
                    laCot2 = laCot * laCot;
                    
                    % psi
                    psi = Elix(i, 1) - e2 * fiSin(i, 1) * fiCos(i, 1) *(1 + c1 * fiSin2(i, 1) + c2 * fiSin4(i, 1) + c3 * fiSin6(i, 1));
                    psiSin = sin(psi);
                    psiCos = cos(psi);
                    psiTan = psiSin / psiCos ;
                    psiSin2 = psiSin * psiSin ;
                    
                    % omega
                    ome = atan(psiTan / laCos);
                    
                    % sigma
                    if laSin ~= 0
                        sigSin =laSin * psiCos;
                        sig = asin(sigSin);
                    else
                        sigSin = 0;
                        sig = 0;
                    end
                    
                    sigSin2 = sigSin * sigSin;
                    sigSin4 = sigSin2 * sigSin2;
                    sigCos2 = 1 - sigSin2;
                    sigCos4 = sigCos2 * sigCos2;
                    
                    % chi
                    chi = sig / 2 + pi / 4;
                    chiTan = tan(chi);
                    chiLog = log(chiTan);
                    
                    % constants
                    aa = psiSin * psiCos * laCos * (1 + sigSin2) / sigCos4;
                    bb = sigSin * (sigCos2 - 2 * psiSin2) / sigCos4;
                    
                    if laCot ~= 0
                        a1 = (psiSin2 - sigSin4 * laCot2)/ sigCos4;
                        b1 = 2 * sigSin2 * psiSin * laCot / sigCos4;
                    else
                        a1 = psiSin2 / sigCos4 ;
                        b1 = 0;
                    end
                    a2 = a1 * a1 - b1 * b1 ;
                    b2 = 2 * a1 * b1 ;
                    a3 = a1 * a2 - b1 * b2 ;
                    b3 = a1 * b2 + b1 * a2 ;
                    rr = k0 - a1 * k1 + a2 * k2 - a3 * k3;
                    tt = b1 * k1 - b2 * k2 + b3 * k3;
                    
                    % X/Y coordinates
                    xx = k * ome + aa * rr + bb * tt;
                    yy = k * chiLog + bb * rr - aa * tt;
                    
                    % North and East
                    if north_south(i, 1) == 1
                        NoEq = 0;
                    else
                        NoEq = 10000000;
                    end
                    
                    north(i, 1) = NoEq + xx;
                    east(i, 1) = EsMc + yy;
                end
            else
                utm_zone = [];
                north = [];
                east = [];
            end
        end
        
        function [loc,  rot_mat] = cart2local(xyz_ref, xyz_baseline)
            % cart2local: from geocentric cartesian baselines (DX) to local coordinates in X0
            %
            % SYNTAX
            %   [loc,  rot_mat] = Coordinates.cart2local(xyz_ref, baseline)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
            loc = (rot_mat * xyz_baseline')';
        end

        function [lat_m, lon_m] = geod2mag(lat_g, lon_g)
            % 'IGRF20120-geomagnetic'
            % from: https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt
            g10 = -29404.8; %1st order IGFR2011 coefficient
            g11 = -1450.9;  %2nd order
            h11 = 4652.5;   %3rd order
            mag.name='IGRF2013-geomagnetic2020.0';
            mag.lam = atan(h11/g11);
            mag.phi    = atan((g11*cos(mag.lam)+h11*sin(mag.lam))/g10);

            % Compose rotation matrix
            Tz=[ cos(mag.lam) sin(mag.lam) 0 ;
                -sin(mag.lam) cos(mag.lam) 0 ;
        	      0       0           1 ];

            Ty=[ cos(mag.phi)  0   -sin(mag.phi) ;
                0       1      0     ;
                sin(mag.phi)  0   cos(mag.phi) ];
            T=Ty*Tz;

            % Degrees to radians, and make into rows.
            lat_g = (pi/2 - lat_g(:)'); % lat to colat
            lon_g = lon_g(:)';

            % Transformation to cartesian coordinates.
            x = sin(lat_g).*cos(lon_g);
            y = sin(lat_g).*sin(lon_g);
            z = cos(lat_g);

            % Rotations of the coordinates.
            tmp = T*[x; y; z];
            xp = tmp(1,:);
            yp = tmp(2,:);
            zp = tmp(3,:);

            % Transformation back to spherical coordinates.
            lat_m = pi/2 - acos(zp./sqrt(xp.^2+yp.^2+zp.^2));
            lon_m = atan2(yp, xp);
        end

        function [lat_m, lon_m] = mag2geod(lat_g, lon_g)
            % 'IGRF20120-geomagnetic'
            % from: https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt
            g10 = -29404.8; %1st order IGFR2011 coefficient
            g11 = -1450.9;  %2nd order
            h11 = 4652.5;   %3rd order
            mag.name='IGRF2013-geomagnetic2020.0';
            mag.lam = atan(h11/g11);
            mag.phi    = atan((g11*cos(mag.lam)+h11*sin(mag.lam))/g10);

            % Compose rotation matrix
            Tz=[ cos(-mag.lam) sin(-mag.lam) 0 ;
                -sin(-mag.lam) cos(-mag.lam) 0 ;
        	      0       0           1 ];

            Ty=[ cos(-mag.phi)  0   -sin(-mag.phi) ;
                0       1      0     ;
                sin(-mag.phi)  0   cos(-mag.phi) ];
            T=Tz*Ty;

            % Degrees to radians, and make into rows.
            lat_g = (pi/2 - lat_g(:)'); % lat to phi (colat)
            lon_g = lon_g(:)';

            % Transformation to cartesian coordinates.
            x = sin(lat_g).*cos(lon_g);
            y = sin(lat_g).*sin(lon_g);
            z = cos(lat_g);

            % Rotations of the coordinates.
            tmp = T*[x; y; z];
            xp = tmp(1,:);
            yp = tmp(2,:);
            zp = tmp(3,:);

            % Transformation back to spherical coordinates.
            lat_m = pi/2 - acos(zp./sqrt(xp.^2+yp.^2+zp.^2));
            lon_m = atan2(yp, xp);
        end
        
        function [xyz, rot_mat] = local2cart(xyz_ref, local_baseline)
            % local2cart: from local coordinates in X0 to geocentric cartesian baselines (DX)
            %
            % SYNTAX
            %   [loc,  rot_mat] = Coordinates.local2cart(xyz_ref, local_baseline)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)]';
            xyz = (rot_mat * local_baseline')';
        end
        
        function [loc,  rot_mat] = cart2loca(xyz_ref, xyz_baseline)
            [loc, rot_mat] = Coordinates.cart2local(xyz_ref, xyz_baseline); % maybe the function was used with the wrong name
        end
        
        function [xyz_baseline, rot_mat] = loca2cart(xyz_ref, loc)
            % loca2cart: from local coordinates in X0 to geocentric cartesian baselines (DX)
            %
            % SYNTAX
            %   [baseline, rot_mat] = Coordinates.loca2cart(xyz_ref, loca)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)]';
            xyz_baseline = (rot_mat * loc)';
        end

        function [lat, lon] = utm2geod(easting, northing, utmzone)
            % Convert UTM coordinates to latitude and longitude
            %
            % Inputs:
            %   easting:    UTM easting in meters
            %   northing:   UTM northing in meters
            %   utmzone
            %       - zone:       UTM zone (integer between 1 and 60)
            %       - hemisphere: UTM hemisphere ('N' or 'S')
            %
            % Outputs:
            %   lat:        Latitude in decimal degrees
            %   lon:        Longitude in decimal degrees
            %
            % Outputs:
            %    Lat: Latitude vector.   Degrees.  +ddd.ddddd  WGS84
            %    Lon: Longitude vector.  Degrees.  +ddd.ddddd  WGS84            
            
            % Argument checking
            %
            narginchk(3, 3); %3 arguments required
            if (size(utmzone,1) == 1)
                utmzone = repmat(utmzone, length(easting),1);
            end
            if (size(utmzone,2) ~= 4)
                error('utmzone should be a vector of strings like "30 T"');
            end

            % Check if utmzone is a vector of strings like "30 T", not "30 t"
            if any(utmzone(:,4) > 'X' | utmzone(:,4) < 'C')
                fprintf('utm2deg: Warning utmzone should be a vector of strings like "30 T", not "30 t"\n');
            end

            % Determine hemisphere
            hemis = repmat('N', size(utmzone,1), 1);
            idx = utmzone(:,4) <= 'M';
            hemis(idx) = 'S';

            % Compute constants
            sa = 6378137.000000 ; sb = 6356752.314245;
            e2 = (((sa ^ 2) - (sb ^ 2)) ^ 0.5) / sb;
            e2cuadrada = e2 ^ 2;
            c = (sa ^ 2) / sb;

            % Compute X and Y values
            X = easting - 500000;
            Y = northing;
            Y(hemis == 'S') = northing(hemis == 'S') - 10000000;

            % Compute latitude and zone values
            S = ((str2double(utmzone(:,1:2)) * 6) - 183);
            lat =  Y ./ (6366197.724 * 0.9996);
            v = (c ./ ((1 + (e2cuadrada .* (cos(lat)) .^ 2))) .^ 0.5) * 0.9996;
            a = X ./ v;

            % Compute intermediate variables
            a1 = sin(2 * lat);
            a2 = a1 .* (cos(lat)) .^ 2;
            j2 = lat + (a1 / 2);
            j4 = ((3 * j2) + a2) / 4;
            j6 = ((5 * j4) + (a2 .* (cos(lat)) .^ 2)) / 3;
            alfa = (3 / 4) * e2cuadrada;
            beta = (5 / 3) * alfa .^ 2;
            gama = (35 / 27) * alfa .^ 3;
            Bm = 0.9996 * c .* (lat - alfa .* j2 + beta .* j4 - gama .* j6);
            b = (Y - Bm) ./ v;
            Epsi = ((e2cuadrada .* a .^ 2) / 2) .* (cos(lat)) .^ 2;
            Eps = a .* (1 - (Epsi / 3));
            nab = (b .* (1 - Epsi)) + lat;
            senoheps = (exp(Eps) - exp(-Eps)) ./ 2;
            Delt = atan(senoheps ./ (cos(nab)));
            TaO = atan(cos(Delt) .* tan(nab));
            lon = (Delt *(180 / pi)) + S;
            lat = (lat + (1 + e2cuadrada .* (cos(lat).^ 2) - (3 / 2) .* e2cuadrada .* sin(lat) .* cos(lat) .* (TaO - lat)) .* (TaO - lat)) .* (180 / pi);
        end
        
        function fh = showCompareENU(coo_list)
            % Plot East North Up coordinates
            %
            % SYNTAX
            %   fh = Coordinates.showCoordinatesENU(coo_list);
            
            fh = coo_list.showPositionENU();
        end
        
        function bslCompare(coo_list_0, coo_list_1, id_ref, spline_base, outlier_thr, y_lim)
            % Compare two sets of receivers coordinates (e.g. coming from two different execution)
            %
            % INPUT
            %   coo_list_0           first set of GNSS receivers or coordinates
            %   coo_list_1           second set of GNSS receivers or coordinates
            %   id_ref               id of the reference receiver
            %                        DEFAULT: last receiver
            %   spline_base          spline base in days (for removing splines)
            %                        DEFAULT: 365/4
            %   outlier_threshold    threshold on the baseline difference for the
            %                        spline and STD computation [mm]
            %                        DEFAULT: 2 mm
            %   y_lim                3 x 2 ylimits
            %
            % SINTAX
            %   Coordinates.bslCompare(rec_list_0, rec_list_1, id_ref, spline_base, outlier_thr)
            %
            % EXAMPLE
            %   Coordinates.bslCompare(nomp.core.rec, core.rec, 7, 365/4, 2);
            
            if nargin < 3 || isempty(id_ref)
                id_ref = numel(coo_list_0); % set to the last;
            end
            if nargin < 4 || isempty(spline_base)
                spline_base = 365/4;
            end
            % Set outlier treshold (on baseline differrence);
            if nargin < 5 || isempty(outlier_thr)
                outlier_thr = 2;
            end
            
            if isa(coo_list_0, 'Coordinates')
                coo_ref0 = coo_list_0(id_ref);
            else
                coo_ref0 = coo_list_0(id_ref).getPos;
            end
            if isa(coo_list_1, 'Coordinates')
                coo_ref1 = coo_list_1(id_ref);
            else
                coo_ref1 = coo_list_1(id_ref).getPos;
            end
            if nargin < 6 || isempty(y_lim)
                y_lim = [-5 5; -10 10; -10 10];
            elseif spline_base == 0
                spline_base = -1;
            end
            
            Core.getLogger.addMonoMessage('\nProcessing gain - solution 0 vs 1\n----------------------------------------');
            
            for r = setdiff(1 : numel(coo_list_1), id_ref)
                if isa(coo_list_0, 'Coordinates')
                    coo0 = coo_list_0(r);
                else
                    coo0 = coo_list_0(r).getPos;
                end
                if isa(coo_list_1, 'Coordinates')
                    coo1 = coo_list_1(r);
                else
                    coo1 = coo_list_1(r).getPos;
                end
                
                % Extract the baselines
                [t_comm, idx1, idx2] = intersect(round(coo_ref0.time.getRefTime(coo0.time.first.getMatlabTime)),round(coo0.time.getRefTime(coo0.time.first.getMatlabTime)));
                t0 = coo0.time.first.getMatlabTime + t_comm/86400;
                enu_diff0 = Coordinates.cart2local(median(coo_ref0.xyz,1,'omitnan'),coo0.xyz(idx2,:) - coo_ref0.xyz(idx1,:) )*1e3;
                id_ok = not(any(isnan(enu_diff0), 2)); % discard NaN
                t0 = t0(id_ok);
                enu_diff0 = enu_diff0(id_ok,:);
                
                [t_comm, idx1, idx2] = intersect(round(coo_ref1.time.getRefTime(coo1.time.first.getMatlabTime)),round(coo1.time.getRefTime(coo1.time.first.getMatlabTime)));
                t1 = coo1.time.first.getMatlabTime + t_comm/86400;
                enu_diff1 = Coordinates.cart2local(median(coo_ref1.xyz,1,'omitnan'),coo1.xyz(idx2,:) - coo_ref1.xyz(idx1,:) )*1e3;
                id_ok = not(any(isnan(enu_diff1), 2)); % discard NaN
                t1 = t1(id_ok);
                enu_diff1 = enu_diff1(id_ok,:);
                
                % Sync the baselines
                [t_comm, idx0, idx1] = intersect(round(t0*86400)/86400,round(t1*86400)/86400, 'stable');
                enu_diff = enu_diff0(idx0,:) - enu_diff1(idx1,:);
                enu_diff = bsxfun(@minus, enu_diff, median(enu_diff, 'omitnan'));
                id_ok = abs(enu_diff - median(enu_diff, 'omitnan')) < outlier_thr;
                
                fh = figure; Core_UI.beautifyFig(fh); drawnow
                % Plot the baseline difference ----------------------------------------
                subplot(3,1,1);
                tmp = bsxfun(@minus, enu_diff, [-20 0 20]);
                plotSep(t_comm, tmp, '.-', 'MarkerSize', 8.5, 'LineWidth', 1.6);
                for c = 1 : 3
                    std_enu(:,c) = std(enu_diff(id_ok(:,c), c), 'omitnan');
                end
                legend(sprintf('East (%.2f mm)', std_enu(1)), sprintf('North (%.2f mm)', std_enu(2)), sprintf('Up (%.2f mm)', std_enu(3)), 'location', 'EastOutside');
                ylim(y_lim(1,:));
                xlim(minMax(t_comm));
                setTimeTicks();
                ax(1) = gca;
                ylabel(sprintf('Baseline\ndifference'));
                grid minor
                if isempty(coo1.name)
                    trg_rec_name = sprintf('%d', r);
                else
                    trg_rec_name = coo1.name;
                end
                if isempty(coo_ref1.name)
                    ref_rec_name = sprintf('%d', id_ref);
                else
                    ref_rec_name = coo_ref1.name;
                end
                title(sprintf('Baseline %s - %s\\fontsize{5} \n', strrep(trg_rec_name,'_','\_'), strrep(ref_rec_name,'_','\_')), 'FontSize', 16);
                
                
                
                % Compute common reduction by spline
                tmp0 = enu_diff0(idx0,:) - median(enu_diff0(idx0,:), 'omitnan');
                tmp1 = enu_diff1(idx1,:) - median(enu_diff1(idx1,:), 'omitnan');
                splined = nan(size(tmp0,1),3);
                lid_ko0 = false(size(tmp0,1),3);
                lid_ko1 = false(size(tmp0,1),3);
                if spline_base > 0
                    for c = 1 : 3
                        ttmp = t_comm(id_ok(:,c));
                        [filtered0, lid_ko0(id_ok(:,c),c), trend0, splined0] = Coordinates.cooFilter([ttmp, tmp0(id_ok(:,c),c)], 0.8, 7,[28, 28, 28]);
                        [filtered1, lid_ko1(id_ok(:,c),c), trend1, splined1] = Coordinates.cooFilter([ttmp, tmp1(id_ok(:,c),c)], 0.8, 7,[28, 28, 28]);
                        splined(id_ok(:,c),c) = mean([splined0 splined1],2);
                    end
                else
                    splined = nan2zero(splined);
                end
                id_ok = not(lid_ko0 | lid_ko1);
                
                % Plot the baseline (filtered by spline) of the solution with no MP ---
                subplot(3,1,2);
                tmp = enu_diff0(idx0,:) - median(enu_diff0(idx0,:), 'omitnan');
                tmp0 = enu_diff0 - median(enu_diff0(idx0,:), 'omitnan');
                splined0 = zeros(size(tmp0));
                if spline_base > 0
                    for c = 1 : 3
                        ttmp = t_comm(id_ok(:,c));
                        [~, ~, ~, splined0(:,c)] =       splinerMat(ttmp, splined(id_ok(:,c),c), spline_base, 1e-8, t0);
                    end
                end
                
                tmp = tmp - splined; % remove splines
                tmp0 = tmp0 - splined0; % remove splines
                tmp0 = bsxfun(@minus, tmp0, [-20 0 20]);
                plotSep(t0, tmp0, '.-', 'MarkerSize', 8.5, 'LineWidth', 1.6);
                
                for c = 1 : 3
                    std_enu0(:,c) = std(tmp(id_ok(:,c), c), 'omitnan');
                end
                %plotSep(t_comm, splined, '.-', 'MarkerSize', 1, 'LineWidth', 1, 'Color', 'k');
                legend(sprintf('East (%.2f mm)', std_enu0(1)), sprintf('North (%.2f mm)', std_enu0(2)), sprintf('Up (%.2f mm)', std_enu0(3)), 'location', 'EastOutside');
                if spline_base ~= 0
                    ylim(y_lim(2,:));
                end
                xlim(minMax(t_comm));
                setTimeTicks();
                ax(2) = gca;
                ylabel(sprintf('Baseline 0\n(reduced)'));
                grid minor
                
                % Plot the baseline (filtered by spline) of the solution with MP ------
                subplot(3,1,3);
                tmp = enu_diff1(idx1,:) - median(enu_diff1(idx1,:), 'omitnan');
                tmp1 = enu_diff1 - median(enu_diff1(idx1,:), 'omitnan');
                splined1 = zeros(size(tmp1));
                if spline_base > 0
                    for c = 1 : 3
                        ttmp = t_comm(id_ok(:,c));
                        [~, ~, ~, splined1(:,c)] =       splinerMat(ttmp, splined(id_ok(:,c),c), spline_base, 1e-8, t1);
                    end
                end
                tmp = tmp - splined; % remove splines
                tmp1 = tmp1 - splined1; % remove splines
                tmp1 = bsxfun(@minus, tmp1, [-20 0 20]);
                plotSep(t1, tmp1, '.-', 'MarkerSize', 8.5, 'LineWidth', 1.6);
                
                for c = 1 : 3
                    std_enu1(:,c) = std(tmp(id_ok(:,c), c), 'omitnan');
                end
                legend(sprintf('East (%.2f mm)', std_enu1(1)), sprintf('North (%.2f mm)', std_enu1(2)), sprintf('Up (%.2f mm)', std_enu1(3)), 'location', 'EastOutside');
                if spline_base ~= 0
                    ylim(y_lim(3,:));
                end
                xlim(minMax(t_comm));
                setTimeTicks();
                ax(3) = gca;
                ylabel(sprintf('Baseline 1\n(reduced)'));
                grid minor
                
                Core_UI.beautifyFig(fh);
                Core_UI.addBeautifyMenu(fh);
                Core_UI.addExportMenu(fh);
                Core_UI.addLineMenu(fh)
                linkaxes(ax, 'x');
                
                Core.getLogger.addMonoMessage(sprintf('Baseline %d - %d) %5.2f %% %5.2f %% %5.2f %%', r, id_ref, (100*((std_enu0 - std_enu1) ./ std_enu0))));
            end
        end
        
        function [data, lid_ko, running_mean, spline, jump_list] = cooFilter(data, robustness_perc, n_sigma, spline_base)
            % Returns the data removing outliers (spikes)
            %
            % INPUT:
            %   data                column array of values
            %   robustness_perc     maximum percentage of date with no outliers
            %
            % SYNTAX:
            %   [data, id_ko] = Coordinates.cooFilter(data, robustness_perc)
            
            if any(data(:,end))
                n_data = size(data,1);
                if nargin < 2
                    robustness_perc = 0.8;
                end
                if nargin < 3
                    n_sigma = 6;
                end
                if nargin < 4 || numel(spline_base) ~= 3
                    spline_base = [28, 21, 7]*86400;
                    if size(data,2) >= 2
                        time = (data(:,1) - data(1))*86400;
                        rate = round(median(diff(time)));
                        if rate > 43200
                            spline_base(3) = max(21, spline_base(3));
                        end
                    end
                end
                flag_time = false;
                idf = [];
                if size(data,2) >= 2
                    flag_var = true;
                    if size(data,2) >= 3
                        data_var = data(:,3);
                    else
                        flag_var = false;
                    end
                    
                    % Suppose regularly sampled data, fill missing epochs with nan
                    time = (data(:,1) - data(1))*86400;
                    rate = round(median(diff(time)));
                    if isnan(rate)
                        rate = 8640.0;
                    end
                    time_full = linspace(time(1), time(end), round((time(end) - time(1)) / rate + 1))';
                    [~, idf, idr] = intersect(round((time_full - rate/2)/rate, 3), round((time - rate/2)/rate, 3));
                    tmp = data(:,2);
                    data = nan(numel(time_full), 1);
                    if numel(idr) < numel(tmp)
                        Core.getLogger.addWarning('StrongFilter is loosing some observations out of sync');
                    end
                    data(idf) = tmp(idr);
                    flag_time = 1;
                else
                    rate = 86400;
                    flag_var = false;
                    idf = 1:numel(data);
                    idr = idf;
                end
                
                % Compute a trend "robust" using the robustness_perc of data
                % [tmp, running_mean] = strongDeTrend(data, robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
                
                if 86400/rate > 24
                    rate = max(3600, rate * 4); % for sub hourly solution it is better to use smaller windows
                end
                
                if flag_var
                    [jump_list, lid_ko, tmp, running_mean] = getJumps([data(idf) data_var(idr)], 86400/rate);
                    if sum(lid_ko)/numel(lid_ko) > 0.75
                        jump_list = 0;
                        lid_ko = data_var(idr) > 1e5 | abs(data(idf)) > 1e4;
                        tmp(not(lid_ko)) = data(idf(not(lid_ko)));
                        running_mean = movmedian(tmp, 5, 'omitnan');
                    end
                else
                    [jump_list, lid_ko, tmp, running_mean] = getJumps(data(idf), 86400/rate);
                end
                if any(jump_list)
                    % Dejump tmp
                    dejump = diff(tmp);
                    dejump(noNaN(zero2nan(jump_list))) = 0;
                    tmp = [0; cumsum(dejump)] + tmp(1);
                    % remove jumps from running mean
                    dejump = diff(running_mean); 
                    dejump(noNaN(zero2nan(jump_list))) = 0; 
                    dejump = [0; cumsum(dejump)] + running_mean(1);
                    trend = Core_Utils.interp1LS(idf, dejump, 2);
                    clear dejump;
                    % dejump data
                    data  = data(idf) - trend;
                    dejump = diff(data);
                    jump_reduction = dejump;
                    jump_reduction(setdiff(1:length(jump_reduction), noNaN(zero2nan(jump_list)))) = 0;
                    jump_reduction = [0; cumsum(jump_reduction)];
                    dejump(noNaN(zero2nan(jump_list))) = 0;
                    data = [0; cumsum(dejump)] + data(1);
                else
                    jump_reduction = 0;
                    trend = Core_Utils.interp1LS(idf, running_mean, 2);
                    data  = data(idf) - trend;
                end                
                tmp = tmp - trend;
                
                if any(tmp) && flag_time && (numel(data) > 4)
                    if (numel(tmp) > 11)
                        spline_base = max(1,min(floor(time(end)-time(1)), spline_base)); % minimumum spline => a day
                        warning off;
                        % Perform a bit of outlier detection before computing splines
                        thr = 6 * perc(abs(tmp), 0.8);
                        
                        
                        % Computer long splines (reduce the signal, montly splines)
                        if flag_var
                            % ok values are within thr range with variations less than thr
                            % or within 6*thr with variations less than thr/6
                            lid_ok = not(lid_ko);
                            [~, ~, ~, long_spline] = splinerMat(time(lid_ok), [data(lid_ok) data_var(lid_ok)], spline_base(1), 1e-5, time); % long splines
                        else
                            lid_ok = not(lid_ko);
                            [~, ~, ~, long_spline] = splinerMat(time(lid_ok), [data(lid_ok) tmp(lid_ok).^2], spline_base(1), 1e-5, time); % long splines
                        end
                        
                        % Keep in tmp the reduced value
                        
                        tmp = data - long_spline(idr);
                        
                        % Compute medium splines (reduce the signal weekly splines)
                        [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(lid_ok) abs(tmp(lid_ok))], spline_base(2), 1e-5, time); % medium splines
                        [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(lid_ok) abs(tmp(lid_ok) - spline(lid_ok))], spline_base(2), 1e-5, time); % medium splines
                        
                        % These are the medium long frequencies, I reduce the signal so that the interpolation will be more stable
                        long_spline = long_spline + spline;
                        
                        % Keep in tmp the reduced value
                        tmp = data - long_spline(idr);
                        
                        % Remove high frequencies
                        thr = n_sigma * min(strongStd(tmp, robustness_perc), perc(abs(tmp - median(tmp, 'omitnan')), robustness_perc));
                        lid_ok = abs(tmp) < thr;
                        if sum(lid_ok) > 2
                            if flag_var
                                [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(lid_ok) data_var(lid_ok)], spline_base(3), 1e-2, time); % short splines
                            else
                                [~, ~, ~, spline] = splinerMat(time(lid_ok), [tmp(lid_ok) tmp(lid_ok).^2], spline_base(3), 1e-2, time); % short splines
                            end
                        end
                        warning on;
                        
                        spline = spline(idr) + long_spline(idr);
                        tmp = data - spline;
                        spline = spline + trend;
                    else
                        lid_ok =  true(size(tmp));
                        spline = trend;
                    end
                else
                    lid_ok =  true(size(tmp));
                    spline = trend;
                end
                
                % Outlier detection based on the interpolation
                thr = n_sigma * min(strongStd(tmp, robustness_perc), perc(abs(tmp - median(tmp, 'omitnan')), robustness_perc));
                if flag_var
                    lid_ko = abs(tmp) > thr | data_var(idr) > 2*perc(data_var(lid_ok),0.975);
                else
                    lid_ko = abs(tmp) > thr;
                end
                
                data = data + trend;
                
                %lid_ko = lid_ok | false; % DEBUGG <======================================================
                
                % figure; plot(data, 'Color', [0.5 0.5 0.5]);
                data(lid_ko) = nan;
                
                if n_data > numel(idr)
                    tmp = nan(n_data, 1);
                    tmp(idr) = data;
                    data = tmp;
                    
                    tmp = true(n_data, 1);
                    tmp(idr) = lid_ko;
                    lid_ko = tmp;
                    
                    tmp = nan(n_data, 1);
                    tmp(idr) = running_mean;
                    running_mean = tmp;
                    
                    tmp = nan(n_data, 1);
                    tmp(idr) = spline;
                    spline = tmp;
                end
                data = data + jump_reduction;
                % hold on; plot(data, '.-b', 'LineWidth', 2)
                % plot(tmp,'g');
            else
                % no data
                jump_list = 0;
                data = data(:,end);
                lid_ko = true(size(data));
                running_mean = data;
                spline = data;
            end
        end        
    end
    
    % =========================================================================
    %%    TESTS
    % =========================================================================
    
    methods (Static, Access = 'public')
        function test()
            % Testing function, tests some basic transformations
            % To be done
            %
            % SYNTAX
            %   test()
            log = Core.getLogger();
            
            log.addMessage('Testing Class Coordinates');
            tic;
            pos_diff = 0;
            
            if pos0 == pos1
                log.addStatusOk('Passed');
            else
                log.addWarning(sprintf('Difference greater than 0.2 ms: %e',t_diff));
            end
            toc
        end
    end
    
    % =========================================================================
    % %   UTILITIES
    % =========================================================================
    
    methods (Static, Access = 'public')
        function [el, d, az] = getSigthElevation(phi1, lambda1, r1, phi2, lambda2, r2)
            psi = acosd(cosd(phi1) .* cosd(phi2) .* cosd(lambda1 - lambda2) + sind(phi1) .* sind(phi2));
            d = r1 * tand(psi);
            h = cosd(psi).*r2 - r1;
            el = atan2(h, d) * 180/pi;
            if nargout == 3
                az = atan2(sind(lambda2-lambda1) .* cosd(phi2), cosd(phi1).*sind(phi2) - sind(phi1)*cosd(phi2).*cosd(lambda2-lambda1))*180/pi;
                az = mod(az,360);
            end
        end

        function [d, cospsi] = sphericaldist(phi1, lambda1, r1, phi2, lambda2, r2)
            cospsi = cosd(phi1) .* cosd(phi2) .* cosd(lambda1 - lambda2) + sind(phi1) .* sind(phi2);
            d = (r1.^2 + r2.^2 - 2 .* r1 .* r2 .* cospsi) .^ (1 / 2);
        end
    end

end
