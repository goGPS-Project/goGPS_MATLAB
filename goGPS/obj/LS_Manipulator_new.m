%   CLASS LS_Manipulator_new
% =========================================================================
%
% DESCRIPTION
%   Manipulate least squares system in a new and more rational way w.r.t
%   the old class
%
% EXAMPLE
%   LSM = LS_Manipulator();
%
% SEE ALSO
%   - Least_Square
% FOR A LIST OF CONSTANTs and METHODS use doc Main_Settings

%--------------------------------------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0 beta 4 ION
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Giulio Tagliaferro
%  Contributors:     Andrea Gatti
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
classdef LS_Manipulator_new < handle
    
    % ==================================================================================================================================================
    %% Parameter columns id (order)
    % ==================================================================================================================================================
    
    properties (Constant)
        PAR_REC_X = 1;
        PAR_REC_Y = 2;
        PAR_REC_Z = 3;
        PAR_REC_EB = 4;
        PAR_AMB = 5;
        PAR_REC_CLK = 6;
        PAR_TROPO = 7;
        PAR_TROPO_N = 8;
        PAR_TROPO_E = 9;
        PAR_TROPO_V = 10;
        PAR_SAT_CLK = 11;
        PAR_ANT_MP = 12;
        PAR_IONO = 13;
        PAR_TROPO_S = 14;
        PAR_SAT_X = 15;
        PAR_SAT_Y = 16;
        PAR_SAT_Z = 17;
        PAR_SAT_EB = 18;
        PAR_REC_EB_LIN = 19;
        
        
        CLASS_NAME = {'REC_X', 'REC_Y', 'REC_Z', 'REC_EB', 'AMB', 'REC_CLK', 'TROPO', 'TROPO_N', 'TROPO_E', 'SAT_CLK', 'ANT_MP', 'IONO', 'TROPO_S', 'SAT_X', 'SAT_Y', 'SAT_Z', 'SAT_EB', 'REC_EB_LIN'};
    end
    
    properties
        %%% sparse A , this in more efficient than a sparse(slow to index in matlab) because
        %%% normally all the observation euation have the same numebr of
        %%% entry
        A % observations
        A_idx
        obs
        res
        param_class % class id of the column of A
        time_obs   % epoch of the observation (GPS_Time)
        satellite_obs % goid satellite of the observations
        receiver_obs % receiver of the observations
        azimuth_obs % azimuth of the observations
        elevation_obs % elevation of the observations
        variance_obs % variance of the observations
        obs_codes_id_obs % id of the signal used
        phase_obs % logical to tell if obs are phase or code
        wl_id_obs % id of the wavelength
        outlier_obs % is obs an outlier?
        
        A_pseudo
        A_idx_pseudo
        time_pseudo   % epoch of the pseudo-observation (GPS_Time)
        satellite_pseudo % goid satellite of the pseudo-observation
        receiver_pseudo % receiver of the pseudo-observation
        variance_pseudo % varaince of the pseudo-observation
        
        unique_obs_codes % uniques ids (cell) of the signals / since lot of combinations are possible they will be generated dynamically)
        unique_wl % set of unique wavelength
        rec_xyz % receiver coordinates to be used in
        unique_rec_name % names of the receivers
        unique_sat_goid % unique satellite goids
        cycle_slips = cell(1) % epoch of the cycle slip
        unique_time % unique epoch of the system
        
        time_par   % time of the parameter!!! very important the parameters (within the same class e.g. time for satellite s ) MUST be ordered in chronological order
        param_par  % parametrization of the parameter
        time_min   % ref_time_of the parameter
        rec_par    % receiver of the parameter
        sat_par    % satellite of the parameter
        class_par  % class of the parameter
        obs_codes_id_par  % obs code id fo the parameter
        wl_id_par % wl id of the parameter
        out_par % paramters that are observed only by outlier observation
        phase_par % is pahse coode or both
        
        rec_set % set of receivers
        sat_set % set of satellites
        ch_set  % set of observation codes
        
        N
        idx_rd % idx parameter removed to silve the rank deficency
        ls_parametrization;
        
        x
        
        log
    end
    
    methods
        function this = LS_Manipulator_new()
            this.log = Core.getLogger;
        end
    end
    
    methods
        function addObsEq(this, rec, obs_set, param_selection)
            % add observation equations to the matrices
            %
            % SYNTAX:
            %    this.addObsEq(rec, obs_set)
            n_par = length(param_selection);
            
            if isempty(this.param_class)
                this.param_class = param_selection;
            end
            if nargin < 4
                param_selection = this.param_class;
            end
            
            % --- check all parameters presence and their order -----
            par_rec_x_lid = param_selection == this.PAR_REC_X;
            par_rec_x = sum(par_rec_x_lid) > 0;
            par_rec_y_lid = param_selection == this.PAR_REC_Y;
            par_rec_y = sum(par_rec_y_lid) > 0;
            par_rec_z_lid = param_selection == this.PAR_REC_Z;
            par_rec_z = sum(par_rec_z_lid) > 0;
            
            par_sat_x_lid = param_selection == this.PAR_SAT_X;
            par_sat_x = sum(par_sat_x_lid) > 0;
            par_sat_y_lid = param_selection == this.PAR_SAT_Y;
            par_sat_y = sum(par_sat_y_lid) > 0;
            par_sat_z_lid = param_selection == this.PAR_SAT_Z;
            par_sat_z = sum(par_sat_z_lid) > 0;
            
            par_rec_eb_lid = param_selection == this.PAR_REC_EB;
            par_rec_eb = sum(par_rec_eb_lid) > 0;
            par_rec_eb_lin_lid = param_selection == this.PAR_REC_EB_LIN;
            par_rec_eb_lin = sum(par_rec_eb_lin_lid) > 0;
            par_sat_eb_lid = param_selection == this.PAR_SAT_EB;
            par_sat_eb = sum(par_sat_eb_lid) > 0;
            
            par_amb_lid = param_selection == this.PAR_AMB;
            par_amb = sum(par_amb_lid) > 0;
            
            par_rec_clk_lid = param_selection == this.PAR_REC_CLK;
            par_rec_clk = sum(par_rec_clk_lid) > 0;
            
            par_sat_clk_lid = param_selection == this.PAR_SAT_CLK;
            par_sat_clk = sum(par_sat_clk_lid) > 0;
            
            par_tropo_lid = param_selection == this.PAR_TROPO;
            par_tropo = sum(par_tropo_lid) > 0;
            
            par_tropo_e_lid = param_selection == this.PAR_TROPO_E;
            par_tropo_e = sum(par_tropo_e_lid) > 0;
            
            par_tropo_n_lid = param_selection == this.PAR_TROPO_N;
            par_tropo_n = sum(par_tropo_n_lid) > 0;
            
            par_tropo_v_lid = param_selection == this.PAR_TROPO_V;
            par_tropo_v = sum(par_tropo_v_lid) > 0;
            
            par_tropo_s_lid = param_selection == this.PAR_TROPO_S;
            par_tropo_s = sum(par_tropo_s_lid) > 0;
            
            par_iono_lid = param_selection == this.PAR_IONO;
            par_iono = sum(par_iono_lid) > 0;
            
            % ---- add the recever to the recibers
            if Core_Utils.findAinB(rec.parent.getMarkerName,this.unique_rec_name) == 0
                this.unique_rec_name{end+1} = rec.parent.getMarkerName;
                this.rec_xyz = [this.rec_xyz; rec.getMedianPosXYZ];
                r = size(this.rec_xyz,1);
            else
                r = Core_Utils.findAinB(rec.parent.getMarkerName,this.unique_rec_name);
            end
            
            % --- remove syntetic measurement form the observations -----
            [synt_obs, xs_loc] = rec.getSyntTwin(obs_set);
            xs_loc = zero2nan(xs_loc);
            diff_obs = nan2zero(zero2nan(obs_set.obs) - zero2nan(synt_obs));
            n_obs = sum(sum(diff_obs ~= 0));
            
            % set up requested number of parametrs
            n_epochs = size(obs_set.obs, 1);
            n_stream = size(diff_obs, 2); % number of satellites present in the observation set
            
            
            
            % add signal to the unique
            u_obs_code = unique(cellstr(obs_set.obs_code));
            for i = 1 : length(u_obs_code)
                if Core_Utils.findAinB(u_obs_code(i),this.unique_obs_codes) == 0
                    this.unique_obs_codes{end+1} = u_obs_code{i};
                end
            end
            obs_code_id = Core_Utils.findAinB(cellstr(obs_set.obs_code),this.unique_obs_codes);
            
            
            % add wavelength to the unique
            u_wl = unique(round(obs_set.wl*1e15))/1e15;
            this.unique_wl = unique([this.unique_wl; u_wl]*1e15)/1e15;
            [~, obs_wl_id] = ismember(round(obs_set.wl*1e7), round(this.unique_wl*1e7));
            
            
            mfw_on =  sum(param_selection == this.PAR_TROPO | param_selection == this.PAR_TROPO_E | param_selection == this.PAR_TROPO_N | param_selection == this.PAR_TROPO_V) > 0;
            % get the mapping function for tropo
            if mfw_on
                id_sync_out = obs_set.getTimeIdx(rec.time.first, rec.getRate);
                [~, mfw] = rec.getSlantMF(id_sync_out);
                mfw(mfw  > 60 ) = nan;
                %mfw = mfw(id_sync_out,:); % getting only the desampled values
            end
            
            % check whivh observations are phase ones
            phase_s = obs_set.obs_code(:,2) == 'L'; % first char is the system second MUST be the bervation type (Phase pseudorance doppler snr)
            
            % initliaze the matrices
            A = zeros(n_obs, n_par);
            [obs,satellite_obs, azimuth_obs, elevation_obs, variance_obs, wl_obs] = deal(zeros(n_obs, 1));
            [ obs_codes_id_obs,  wl_obs] = deal(zeros(n_obs, 1,'uint8'));
            [ satellite_obs] = deal(zeros(n_obs, 1,'uint16'));
            
            phase_obs = false(n_obs,1);
            time_obs = GPS_Time();
            obs_count = 1;
            if length(this.cycle_slips) < r || isempty(this.cycle_slips{r})
                this.cycle_slips{r} = {};
                this.cycle_slips{r} = cell(max([this.unique_sat_goid obs_set.go_id']),1);
                for s = 1 : max([this.unique_sat_goid obs_set.go_id'])
                    this.cycle_slips{r}{s} = cell(length(this.unique_obs_codes),1);
                end
            end
            iono_const = GPS_SS.L_VEC(1)^2;
            % fill the A matrix per satellite
            for s = 1 : n_stream
                id_ok_stream = diff_obs(:, s) ~= 0; % check observation existence -> logical array for a "s" stream
                if any(id_ok_stream)
                    %---- some quantities related to the stream
                    obs_stream = diff_obs(id_ok_stream, s);
                    n_obs_stream = length(obs_stream);
                    lines_stream = obs_count + (0:(n_obs_stream - 1));
                    
                    el_stream = obs_set.el(id_ok_stream, s) / 180 * pi;
                    az_stream = obs_set.az(id_ok_stream, s) / 180 * pi;
                    s_go_id   = obs_set.go_id(s);
                    s_s_id    = obs_code_id(s);
                    wl_id     = obs_wl_id(s);
                    if mfw_on
                        mfw_stream = mfw(id_ok_stream, s_go_id);
                    end
                    this.unique_sat_goid = unique([this.unique_sat_goid  s_go_id]);
                    
                    
                    xs_loc_stream = permute(xs_loc(id_ok_stream, s, :), [1, 3, 2]);
                    los_stream = rowNormalize(xs_loc_stream);
                    
                    %--- Fill Observation related vectors------------
                    obs(lines_stream) = obs_stream;
                    satellite_obs(lines_stream) = s_go_id;
                    variance_obs(lines_stream) =  obs_set.sigma(s)^2;
                    azimuth_obs(lines_stream) =  obs_set.sigma(s)^2;
                    elevation_obs(lines_stream) =  obs_set.sigma(s)^2;
                    obs_codes_id_obs(lines_stream) = s_s_id;
                    phase_obs(lines_stream) = phase_s(s);
                    wl_obs(lines_stream) = wl_id;
                    time_obs.addEpoch(obs_set.time.getEpoch(id_ok_stream).getMatlabTime);
                    
                    % ----------- save the cycle slips
                    if phase_s(s)
                        %obs_set.cycle_slip(,s) = true;
                        if any(obs_set.obs(:,s))
                            obs_set.cycle_slip(find(obs_set.obs(:,s) ~= 0, 1, 'first'), s) = true;
                        end
                        this.cycle_slips{r}{s_go_id}{s_s_id} = obs_set.time.getEpoch(find(obs_set.cycle_slip(:, s) & obs_set.obs(:, s) ~= 0));
                    end
                    
                    % ----------- FILL IMAGE MATRIX ------------
                    % ----------- coordinates ------------------
                    if par_rec_x
                        A(lines_stream, par_rec_x_lid) = - los_stream(:,1);
                    end
                    
                    if par_rec_y
                        A(lines_stream, par_rec_y_lid) = - los_stream(:,2);
                    end
                    
                    if par_rec_z
                        A(lines_stream, par_rec_z_lid) = - los_stream(:,3);
                    end
                    
                    if par_sat_x
                        A(lines_stream, par_rec_x_lid) =  los_stream(:,1);
                    end
                    
                    if par_sat_y
                        A(lines_stream, par_rec_y_lid) =  los_stream(:,2);
                    end
                    
                    if par_sat_z
                        A(lines_stream, par_rec_z_lid) =  los_stream(:,3);
                    end
                    % ----------- electronic bias ------------------
                    if par_rec_eb
                        A(lines_stream, par_rec_eb_lid) = 1;
                    end
                    if par_rec_eb_lin
                        A(lines_stream, par_rec_eb_lin_lid) = 1/obs_set.wl(s);
                    end
                    if par_sat_eb
                        A(lines_stream, par_sat_eb_lid) = 1;
                    end
                    % ----------- Abiguity ------------------
                    if par_amb && phase_s(s)
                        A(lines_stream, par_amb_lid) = obs_set.wl(s);
                    end
                    % ----------- Clock ------------------
                    if par_rec_clk
                        A(lines_stream, par_rec_clk_lid) = 1;
                    end
                    if par_sat_clk
                        A(lines_stream, par_sat_clk_lid) = 1;
                    end
                    % ----------- ZTD ------------------
                    if par_tropo
                        A(lines_stream, par_tropo_lid) = mfw_stream;
                    end
                    % ----------- ZTD gradients ------------------
                    if par_tropo_n || par_tropo_e
                        cotan_term = 1 ./ ( sin(el_stream).*tan(el_stream) + 0.0032);
                        if par_tropo_e
                            A(lines_stream, par_tropo_e_lid) = sin(az_stream) .* cotan_term; % east gradient
                        end
                        if par_tropo_n
                            A(lines_stream, par_tropo_n_lid) = cos(az_stream) .* cotan_term; % noth gradient
                        end
                    end
                    if par_tropo_v
                        A(lines_stream, par_tropo_v_lid) = mfw_stream*rec.h_ellips;
                    end
                    % ----------- Ionosphere delay --------------------
                    if par_iono
                        if phase_s(s)
                            A(lines_stream, par_iono_lid) = - obs_set.wl(s)^2/iono_const;
                        else
                            A(lines_stream, par_iono_lid) =   obs_set.wl(s)^2/iono_const;
                        end
                    end
                    obs_count = obs_count + n_obs_stream;
                end
            end
            this.A = [this.A; A];
            this.obs = [this.obs; obs];
            this.satellite_obs = [this.satellite_obs; satellite_obs];
            this.azimuth_obs = [this.azimuth_obs; azimuth_obs];
            this.elevation_obs = [this.elevation_obs; elevation_obs];
            this.obs_codes_id_obs = [this.obs_codes_id_obs; obs_codes_id_obs];
            this.variance_obs = [this.variance_obs; variance_obs];
            this.phase_obs = [this.phase_obs; phase_obs];
            this.wl_id_obs = [this.wl_id_obs; wl_obs];
            
            if isempty(this.time_obs)
                this.time_obs = time_obs;
            else
                this.time_obs.addEpoch(time_obs.getMatlabTime);
            end
            this.receiver_obs = [this.receiver_obs; r*ones(size(phase_obs))];
        end
        
        function bondParamsGenerateIdx(this, ls_parametrization)
            % bond parameters (splines or other models) and generate idx
            %
            % SYNTAX
            %    this.bondParamGenerateIdx(parametrization)
            this.log = Core.getLogger();
            this.ls_parametrization = ls_parametrization;
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            n_obs = size(this.A,1);
            this.A_idx = zeros(size(this.A),'uint32');
            time_min = min(this.time_obs.getMatlabTime);
            this.time_min = this.time_obs.minimum();
            obs_rate = this.time_obs.getRate; % TBD substitute this quantity with the obesravtion minimum rate
            time_obs = round(this.time_obs.getRefTime(time_min)/obs_rate);
            
            cumulative_idx = 0;
            i_col = 1;
            
            % is the observatio  pahse or code
            phpr_unique_obs_codes = zeros(size(this.unique_obs_codes));
            for o = 1 : length(this.unique_obs_codes)
                if this.unique_obs_codes{o}(2) == 'C'
                    phpr_unique_obs_codes(o) = 1;
                elseif this.unique_obs_codes{o}(2) == 'L'
                    phpr_unique_obs_codes(o) = 2;
                end
            end
            % ---- preallocation to speed up
            ch_lid = false(size(this.A,1),1);
            obs_lid = false(size(this.A,1),1);
            
            
            
            this.time_par = zeros(n_obs,2,'uint32');
            this.param_par = zeros(n_obs,4,'uint8'); % time of the parameter
            this.rec_par = zeros(n_obs,1,'int16'); % receiver of the parameters
            this.sat_par = zeros(n_obs,1,'int8');  % receiver of the parameters
            this.class_par = zeros(n_obs,1,'uint8');  % class of the parameter
            this.obs_codes_id_par= zeros(n_obs,1,'int8');  % obs_code id parameters
            this.wl_id_par= zeros(n_obs,1,'uint8');  % obs_code id parameters
            this.phase_par= zeros(n_obs,1,'uint8');  % phse code or both parameter
            
            this.outlier_obs = false(size(this.obs));
            
            ch_set_old = []; % avoid repating expesnive task
            i_p = 1;
            while i_p <= length(this.param_class)
                has_not_expanded = true;
                col_incr = 0;
                p = this.param_class(i_p);
                [parametriz, opt] = ls_parametrization.getParametrization(p);
                % ----------------- defining receiver sets ---------
                if parametriz(2) == ls_parametrization.SING_REC
                    n_rec_set = size(this.rec_xyz,1);
                    rec_set = num2cell(1:n_rec_set);
                elseif parametriz(2) == ls_parametrization.ALL_REC
                    n_rec_set = 1;
                    rec_set = {1:size(this.rec_xyz,1)};
                elseif parametriz(2) == ls_parametrization.MULTI_REC
                    n_rec_set = length(opt.rec_sets);
                    rec_set = opt.rec_sets;
                end
                % ----------------- defining satellites sets ---------
                if parametriz(3) == ls_parametrization.SING_SAT
                    n_sat_set = length(this.unique_sat_goid);
                    sat_set = num2cell(this.unique_sat_goid);
                elseif parametriz(3) == ls_parametrization.ALL_SAT
                    n_sat_set = 1;
                    sat_set = {this.unique_sat_goid};
                elseif parametriz(3) == ls_parametrization.MULTI_SAT
                    n_sat_set = length(opt.sat_sets);
                    sat_set = opt.sat_sets;
                end
                
                % ----------------- defining channel sets ---------
                % construct channel to frequency mapping
                sig2wl = zeros(size(this.unique_obs_codes));
                sig_p_id = 1:numel(this.unique_obs_codes);
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        sig2wl(c) = this.wl_id_obs(idx1);
                    end
                end
                % construct channel to phase mapping
                sig2phase = zeros(size(this.unique_obs_codes));
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        sig2phase(c) = this.phase_obs(idx1);
                    end
                end
                % construct channel to constellation mapping
                sig2const = char(zeros(size(this.unique_obs_codes)));
                sig2const_num = zeros(size(this.unique_obs_codes));
                
                for c = 1 : length(this.unique_obs_codes)
                    idx1 = find(this.obs_codes_id_obs == c,1,'first');
                    if ~isempty(idx1)
                        ant_id = Core.getConstellationCollector.getAntennaId(this.satellite_obs(idx1));
                        sig2const(c) = ant_id(1);
                        sig2const_num(c) = find(Core.getConstellationCollector.SYS_C == ant_id(1));
                    end
                end
                if parametriz(4) == ls_parametrization.SING_TRACK
                    n_ch_set = numel(this.unique_obs_codes);
                    ch_set = num2cell(uint8(1:n_ch_set));
                elseif parametriz(4) == ls_parametrization.ALL_FREQ
                    n_ch_set = 1;
                    ch_set = {uint8(1:numel(this.unique_obs_codes))};
                elseif parametriz(4) == ls_parametrization.SING_FREQ
                    n_ch_set = length(this.unique_wl);
                    ch_set = {};
                    for c = 1 : n_ch_set
                        ch_set{c} = uint8(sig_p_id(sig2wl == c));
                    end
                elseif parametriz(4) == ls_parametrization.FREQ_CONST
                    u_wl_const = unique(sig2wl*100+sig2const_num)
                    n_ch_set = length(u_wl_const);
                    ch_set = {};
                    for c = 1 : n_ch_set
                        const = rem(u_wl_const(c),100);
                        wl = floor(u_wl_const(c)/100);
                        ch_set{c} = uint8(sig_p_id(sig2const_num == const & sig2wl == wl));
                    end
                elseif parametriz(4) == ls_parametrization.SING_FREQ_BIN
                    n_ch_set = length(this.unique_wl);
                    ch_set = {};
                    % waveleght are binned into theese bands
                    fr_bin = [1176.450*1e6+[-1e5  +1e5];  %  L5a   -> G5 , I5 , J5, E5 , S5
                        1191.795*1e6+[-1e5  +1e5];  %  L5ab  -> E8
                        1202.025*1e6+[-1e5  +1e5];  %  L5c   -> R3
                        1207.140*1e6+[-1e5  +1e5];  %  L5B   -> C7 , E7,
                        1227.600*1e6+[-1e5  +1e5];  %  L2a   -> G2 , J2
                        1246.000*1e6+[-7*7/16*1e6-1e5  +12*7/16*1e6+1e5];  %  L2b   -> R2
                        1268.520*1e6+[-1e5  +1e5];  %  L6a   -> C3
                        1278.750*1e6+[-1e5  +1e5];  %  L6b   -> J6 , E6
                        1561.098*1e6+[-1e5  +1e5];  %  L1a   -> C2
                        1575.42*1e6+[-1e5  +1e5];  %  L1b   -> G1, E1, J1, S1
                        1602.000*1e6+[-7*9/16*1e6-1e5  +12*9/16*1e6+1e5];  %  L1c   -> R1
                        2492.028*1e6+[-1e5  +1e5];  %  L9   -> I9
                        ];
                    wl_bin = fliplr(Core_Utils.V_LIGHT./fr_bin);
                    wl2bnd = zeros(size(this.unique_wl));
                    wl2chid = zeros(size(this.unique_wl));
                    
                    for w = 1 : length(this.unique_wl)
                        b = 1;
                        while wl2bnd(w) == 0
                            if this.unique_wl(w) > wl_bin(b,1) && this.unique_wl(w) < wl_bin(b,2)
                                wl2bnd(w) = b;
                            else
                                b = b+1;
                            end
                        end
                        wl2chid(w) = this.obs_codes_id_obs(find(this.wl_id_obs == w,1,'first'));
                    end
                    u_bnd = unique(wl2bnd);
                    n_ch_set = length(u_bnd);
                    for c = 1 : n_ch_set
                        chids = wl2chid(wl2bnd == u_bnd(c));
                        ch_set{c} = chids;
                    end
                elseif parametriz(4) == ls_parametrization.RULE
                    % ---- evaluate the rule
                    n_ch_set = 0;
                    ch_set = {};
                    for rule = opt.rule
                        cond_sig = true(size(sig_p_id));
                        parts = strsplit(rule{1},':');
                        cond = parts{1};
                        optt = parts{2};
                        conds = strsplit(cond,'&');
                        for c = conds
                            if strcmpi('PSRANGE',c)
                                cond_sig = cond_sig & ~sig2phase;
                            elseif strcmpi('PHASE',c)
                                cond_sig = cond_sig & sig2phase;
                            elseif strcmpi('GLONASS',c)
                                cond_sig = cond_sig & sig2const == 'R';
                            elseif strcmpi('NOT*GLONASS',c)
                                cond_sig = cond_sig & sig2const ~= 'R';
                            else
                                cond_sig(:) = false;
                                this.log.addWarning(sprintf('Unknown option %s in parametrization rule, ignoring parameter',c));
                            end
                        end
                        if any(cond_sig)
                            if str2num(optt) == ls_parametrization.ALL_FREQ
                                n_ch_set = n_ch_set +1;
                                ch_set{n_ch_set} = uint8(sig_p_id(cond_sig));
                            elseif str2num(optt) == ls_parametrization.SING_FREQ
                                u_wl_tmp = unique(sig2wl(cond_sig));
                                n_wl_tmp = length(u_wl_tmp);
                                for c = 1 : n_wl_tmp
                                    n_ch_set = n_ch_set +1;
                                    ch_set{n_ch_set} = uint8(sig_p_id(sig2wl == u_wl_tmp(c) & cond_sig));
                                end
                            elseif str2num(optt) == ls_parametrization.SING_TRACK
                                n_ch_set = n_ch_set +length(sig_p_id(cond_sig));
                                ch_set = [ch_set num2cell(uint8(sig_p_id(cond_sig)))];
                            elseif str2num(optt) == ls_parametrization.FREQ_CONST
                                u_wl_const = unique(sig2wl(cond_sig)*100+sig2const_num(cond_sig));
                                for c = 1 : length(u_wl_const)
                                    n_ch_set = n_ch_set +1;
                                    const = rem(u_wl_const(c),100);
                                    wl = floor(u_wl_const(c)/100);
                                    ch_set{n_ch_set} = uint8(sig_p_id(sig2const_num == const & sig2wl == wl & cond_sig));
                                end
                            end
                            
                        end
                        
                    end
                end
                
                % ------- Now constructing the index, the epoch dependednce will be dealt inside the loop -------------
                for r = 1 : n_rec_set
                    rec_lid = false(size(this.A,1),1);
                    for rr = rec_set{r}
                        rec_lid = rec_lid | this.receiver_obs == rr;
                    end
                    % find an id for the receiver set to keep track of the
                    % parameters if recievr is sigle is simply the receiver
                    % progessince number, otherwise they are negative
                    % number with the index in the receiver set
                    if length(rec_set{r}) == 1
                        r_id = rec_set{r};
                    else
                        r_id = -Core_Utils.findAinB(rec_set{r},this.rec_set);
                        if r_id == 0
                            r_id = -length(this.rec_set) -1;
                            this.rec_set{end+1} = rec_set{r};
                        end
                        
                    end
                    for s = 1 : n_sat_set
                        %                         sat_lid = false(size(this.A,1),1);
                        %                         for ss = sat_set{s}
                        %                             sat_lid = sat_lid | this.satellite_obs == ss;
                        %                         end
                        sat_lid = ismember(this.satellite_obs,sat_set{s});
                        ss = sat_set{s};
                        if length(sat_set{s}) == 1
                            s_id = sat_set{s};
                        else
                            s_id = -Core_Utils.findAinB(sat_set{s},this.sat_set);
                            if s_id == 0
                                s_id = -length(this.sat_set) -1;
                                this.sat_set{end+1} = sat_set{s};
                            end
                            
                        end
                        for f = 1 : n_ch_set
                            if p ~= this.PAR_AMB || this.unique_obs_codes{ch_set{f}}(2) == 'L'
                                %                                 ch_lid(:) = false;
                                %                                 for cc = ch_set{f}
                                %                                     ch_lid = ch_lid | this.obs_codes_id_obs == cc;
                                %                                 end
                                if ~isequal(ch_set_old,ch_set{f})
                                    %ch_lid = ismemberBuiltinTypes(this.obs_codes_id_obs,ch_set{f});
                                    ch_lid = ismember(this.obs_codes_id_obs, ch_set{f});
                                    ch_set_old = ch_set{f};
                                end
                                cc = ch_set{f};
                                if length(ch_set{f}) == 1
                                    ch_id = ch_set{f};
                                else
                                    ch_id = -Core_Utils.findAinB(ch_set{f},this.ch_set);
                                    if ch_id == 0
                                        ch_id = -length(this.ch_set) -1;
                                        this.ch_set{end+1} = ch_set{f};
                                    end
                                end
                                % is a phase pr or both paramter
                                phase_code = 0;
                                u_phpr = unique(phpr_unique_obs_codes(cc));
                                if length(u_phpr) == 1
                                    phase_code = u_phpr;
                                end
                                %---- final observation id
                                obs_lid = rec_lid & sat_lid & ch_lid;
                                if any(obs_lid)
                                    % --- now dealing with epoch dependence ----
                                    cols_tmp = 0;
                                    if parametriz(1) == ls_parametrization.CONST
                                        ep_pgr_id = 1;
                                        n_prg_id = 1;
                                        time_par_tmp = [this.time_obs.getEpoch(obs_lid).minimum.getRefTime(this.time_min.getMatlabTime)  this.time_obs.getEpoch(obs_lid).maximum.getRefTime(this.time_min.getMatlabTime)];
                                    elseif parametriz(1) == ls_parametrization.EP_WISE
                                        ep_id = time_obs(obs_lid);
                                        u_e_tmp = unique(ep_id);
                                        [~,ep_pgr_id] = ismember(ep_id,u_e_tmp);
                                        n_prg_id = length(u_e_tmp);
                                        time_par_tmp = [u_e_tmp*obs_rate zeros(size(u_e_tmp))];
                                    elseif parametriz(1) == ls_parametrization.STEP_CONST
                                        ep_id = time_obs(obs_lid);
                                        ep_pgr_id = zeros(size(ep_id));
                                        n_prg_id = 0;
                                        time_par_tmp = [];
                                        if p == this.PAR_AMB %% ambiguity case is really unique and make sense to treat it separatly
                                            p_s = 1;
                                            if length(this.cycle_slips{rr}{ss}) >= cc
                                                cs = this.cycle_slips{rr}{ss}{cc};
                                                if ~isempty(cs)
                                                    steps = round(cs.getRefTime(time_min)/obs_rate);
                                                    time_par_tmp = zeros(length(steps),2);
                                                    for st = steps'
                                                        lid_maj = ep_id >= st;
                                                        ep_pgr_id(lid_maj) = p_s;
                                                        time_par_tmp(p_s,:) = [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]; %start of the arc
                                                        if p_s > 1
                                                            last_ep = ep_id(find(~lid_maj,1,'last'));
                                                            if ~isempty(last_ep)
                                                                time_par_tmp(p_s-1,2) = last_ep*obs_rate; % end of the arc
                                                            else
                                                                time_par_tmp(p_s-1,2) = ep_id(end)*obs_rate; % end of the arc
                                                            end
                                                        end
                                                        p_s = p_s +1;
                                                    end
                                                    n_prg_id = p_s-1;
                                                else
                                                    ep_pgr_id = 1;
                                                    n_prg_id = 1;
                                                    time_par_tmp = [ep_id(1)*obs_rate ep_id(end)*obs_rate];
                                                    
                                                end
                                            end
                                        else
                                            steps_set = opt.steps_set;
                                            
                                            if parametriz(2) == ls_parametrization.ALL_REC % you can use differents step for step wise satellite dependent paraters
                                                steps = round(steps_set{ss}.getRefTime(time_min)/obs_rate);
                                                p_s = 1;
                                                for st = steps'
                                                    lid_maj = ep_id >= st;
                                                    ep_pgr_id(lid_maj) = p_s;
                                                    time_par_tmp = [time_par_tmp; [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]]; %start of the arc
                                                    if p_s > 1
                                                        time_par_tmp(p_s-1,2) = ep_id(find(~lid_maj,1,'last'))*obs_rate; % end of the arc
                                                    end
                                                    p_s = p_s +1;
                                                end
                                            elseif parametriz(3) == ls_parametrization.ALL_SAT  % you can use differents step for step wise receiver dependent paraters
                                                steps = round(steps_set{rr}.getRefTime(time_min)/obs_rate);
                                                p_s = 1;
                                                for st = steps'
                                                    lid_maj = ep_id >= st;
                                                    ep_pgr_id(lid_maj) = p_s;
                                                    time_par_tmp = [time_par_tmp; [ep_id(find(lid_maj,1,'first'))*obs_rate ep_id(end)*obs_rate]]; %start of the arc
                                                    if p_s > 1
                                                        time_par_tmp(p_s-1,2) = ep_id(find(~lid_maj,1,'last'))*obs_rate; % end of the arc
                                                    end
                                                    p_s = p_s +1;
                                                end
                                            end
                                            n_prg_id = p_s - 1;
                                            
                                        end
                                    elseif parametriz(1) == ls_parametrization.SPLINE_ZERO
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        u_e_tmp = unique(ep_id);
                                        [~,ep_pgr_id] = ismember(ep_id,u_e_tmp);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        n_prg_id = length(u_e_tmp);
                                    elseif parametriz(1) == ls_parametrization.SPLINE_LIN
                                        % ---- colum will be doubled -----
                                        cols_tmp = [ 0 1];
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate)/opt.spline_rate,1);
                                        u_e_tmp = unique([ep_id ep_id+1]);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                        for i_o = cols_tmp;
                                            [~,ep_pgr_id(i_o+1)] = ismember(ep_id+i_o,u_e_tmp);
                                        end
                                        % ----- expand colum of the A matrix
                                        a_col = this.A(obs_lid,i_col);
                                        if has_not_expanded
                                            this.A = [this.A(:,1:(i_col)) zeros(n_obs,1) this.A(:,(i_col+1):end)];
                                            this.A_idx = [this.A_idx(:,1:(i_col)) zeros(n_obs,1,'uint32') this.A_idx(:,(i_col+1):end)];
                                            this.param_class = [this.param_class(1:(i_col-1)); p*ones(2,1); this.param_class((i_col+1):end)];
                                            i_p = i_p + 1;
                                            has_not_expanded = false;
                                        end
                                        this.A = [this.A(:,1:(i_col-1)) zeros(n_obs,2) this.A(:,(i_col+1):end)];
                                        this.A(obs_lid,i_col + cols_tmp) = [a_col.*spline_v(:,1) a_col.*spline_v(:,2)];
                                        n_prg_id = length(u_e_tmp);
                                    elseif parametriz(1) == ls_parametrization.SPLINE_CUB
                                        % ---- colum will be quadrupled -----
                                        cols_tmp = [ 0 1 2 3];
                                        ep_id = floor(time_obs(obs_lid)*obs_rate/opt.spline_rate);
                                        spline_v = Core_Utils.spline(rem(time_obs(obs_lid)*obs_rate,opt.spline_rate)/opt.spline_rate,3);
                                        u_e_tmp = unique([ep_id ep_id+1 ep_id+2 ep_id+3]);
                                        time_par_tmp = [u_e_tmp*opt.spline_rate  (u_e_tmp+1)*opt.spline_rate];
                                        ep_pgr_id = zeros(sum(obs_lid),length(cols_tmp));
                                        for i_o = cols_tmp;
                                            [~,ep_pgr_id(:,i_o+1)] = ismember(ep_id+i_o,u_e_tmp);
                                        end
                                        % ----- expand colum of the A matrix
                                        
                                        if has_not_expanded
                                            this.A = [this.A(:,1:(i_col)) zeros(n_obs,3) this.A(:,(i_col+1):end)];
                                            this.A_idx = [this.A_idx(:,1:(i_col)) zeros(n_obs,3,'uint32') this.A_idx(:,(i_col+1):end)];
                                            this.param_class = [this.param_class(1:(i_col-1)); p*ones(4,1); this.param_class((i_col+1):end)];
                                            i_p = i_p + 3;
                                            has_not_expanded = false;
                                        end
                                        a_col = this.A(obs_lid,i_col);
                                        this.A(obs_lid,i_col + cols_tmp) = [a_col.*spline_v(:,1) a_col.*spline_v(:,2) a_col.*spline_v(:,3) a_col.*spline_v(:,4)];
                                        n_prg_id = length(u_e_tmp);
                                    end
                                    this.A_idx(obs_lid, i_col + cols_tmp) = cumulative_idx + ep_pgr_id;
                                    [u_new_par] = unique(cumulative_idx + ep_pgr_id(:));
                                    
                                    this.time_par(cumulative_idx+(1:n_prg_id),:) =  uint32(time_par_tmp);% time of the parameter
                                    this.param_par(cumulative_idx+(1:n_prg_id),:) = repmat(uint8(parametriz),n_prg_id,1);% time of the parameter
                                    this.rec_par(cumulative_idx+(1:n_prg_id)) = r_id*ones(n_prg_id,1,'int8');  % receiver of the parameters
                                    this.sat_par(cumulative_idx+(1:n_prg_id)) = s_id*ones(n_prg_id,1,'int8');  % receiver of the parameters
                                    this.class_par(cumulative_idx+(1:n_prg_id)) =  p*ones(n_prg_id,1,'uint8');  % class of the parameter
                                    this.obs_codes_id_par(cumulative_idx+(1:n_prg_id)) = int8(ch_id)*ones(n_prg_id,1,'int8');  % obs_code id parameters
                                    this.phase_par(cumulative_idx+(1:n_prg_id)) = uint8(phase_code)*ones(n_prg_id,1,'uint8');  % phase code or both id parameters
                                    
                                    if ch_id > 0
                                        this.wl_id_par(cumulative_idx+(1:n_prg_id)) = uint8(sig2wl(ch_id))*ones(n_prg_id,1,'uint8');  % wavelength id parameters
                                    else
                                        u_wl = unique(sig2wl(this.ch_set{-ch_id})); % check if at least they belong to the same frequency
                                        if length(u_wl) == 1
                                            this.wl_id_par(cumulative_idx+(1:n_prg_id)) = uint8(u_wl)*ones(n_prg_id,1,'uint8');  % wavelength id parameters
                                        else
                                            this.wl_id_par(cumulative_idx+(1:n_prg_id)) = zeros(n_prg_id,1,'uint8');  % wavelength id parameters
                                        end
                                    end
                                    cumulative_idx = cumulative_idx + n_prg_id;
                                    col_incr = max(col_incr,length(cols_tmp));
                                end
                            end
                        end
                    end
                end
                i_col = i_col + col_incr;
                i_p = i_p + 1;
            end
            % free allocated space in excess
            this.time_par((cumulative_idx+1):end,:) =  [];%
            this.param_par((cumulative_idx+1):end,:) =  [];%
            this.rec_par((cumulative_idx+1):end) =  [];%
            this.sat_par((cumulative_idx+1):end) =  [];%
            this.class_par((cumulative_idx+1):end) =  [];%
            this.obs_codes_id_par((cumulative_idx+1):end) =  [];%
            this.wl_id_par((cumulative_idx+1):end) =  [];%
            this.phase_par((cumulative_idx+1):end) =  [];%
            this.out_par =  false(size(this.wl_id_par));%
            
            
            
            
        end
        
        function generateOutliedParIdx(this)
            % considering outlier generate an index of outlied paramters (i.e paramters that are observed only by outlier observation)
            %
            % SYNTAX
            %   this.generateOutliedParIdx();
            if isempty(this.outlier_obs)
                this.outlier_obs = false(size(this.obs));
            end
            out_par = unique(serialize(this.A_idx(this.outlier_obs > 0,:)));
            not_out_par = unique([serialize(this.A_idx(this.outlier_obs == 0,:)); this.A_idx_pseudo(:)]);
            this.out_par = Core_Utils.ordinal2logical(setdiff(out_par, not_out_par),max(max(this.A_idx)));
            
        end
        
        function markShortArcs(this, arc_length)
            % mark arc equal or shorter than arc_length as outlier
            %
            % SYNTAX:
            %    this.markShortArcs(arc_length)
            amb = double(this.A_idx(~this.outlier_obs,this.param_class == this.PAR_AMB));
            [occ,u_amb]=hist(amb,unique(amb));
            to_mark = find(occ <= arc_length);
            for m = to_mark
                idx_o = amb == u_amb(m);
                this.outlier_obs(idx_o) = true;
            end
            
        end
        
        function markSingledObs(this)
            % mark single satellite or receiver epoch as outlier
            %
            % SYNTAX:
            %    this.markSingledObs()
            for s = this.unique_sat_goid
                idx_par = find(this.sat_par == s & this.class_par == this.PAR_SAT_CLK & ~this.out_par);
                idx_obs = this.satellite_obs == s & ~this.outlier_obs;
                num_rec = zeros(length(idx_par),1);
                for r = 1 : length(this.unique_rec_name)
                    [~,idx_par_rec] = ismember(unique(this.A_idx(this.receiver_obs == r &  idx_obs, this.param_class == this.PAR_SAT_CLK)),idx_par);
                    num_rec(noZero(idx_par_rec)) = num_rec(noZero(idx_par_rec)) + 1;
                end
                idx_el = find(num_rec < 2);
                for i = idx_el'
                    obs_par = this.A_idx(:,this.param_class == this.PAR_SAT_CLK) == idx_par(i);
                    this.outlier_obs(obs_par) = true;
                end
            end
            for r = 1 : length(this.unique_rec_name)
                idx_par = find(this.rec_par == r & this.class_par == this.PAR_REC_CLK & ~this.out_par);
                idx_obs = this.receiver_obs == r & ~this.outlier_obs;
                num_sat = zeros(length(idx_par),1);
                for s = this.unique_sat_goid
                    [~,idx_par_rec] = ismember(unique(this.A_idx(this.satellite_obs == s &  idx_obs, this.param_class == this.PAR_REC_CLK)), idx_par);
                    num_sat(noZero(idx_par_rec)) = num_sat(noZero(idx_par_rec)) + 1;
                end
                idx_el = find(num_sat < 2);
                for i = idx_el'
                    obs_par = this.A_idx(:,this.param_class == this.PAR_REC_CLK) == idx_par(i);
                    this.outlier_obs(obs_par) = true;
                end
            end
            
        end
        
        function removeFullRankDeficency(this)
            % solve full rank deficency removing parameters from the
            % estimation
            %
            % SYNTAX:
            %    this.removeFullRankDeficency()
            
            % TODO make a log of what is taken off to solve the rank
            % deficency
            
            % remove two bias per receiver
            n_rec = size(this.rec_xyz,1);
            n_sat = length(this.unique_sat_goid);
            this.idx_rd = []; %empty previous par choosen to solve the rank deficency
            idx_rm = [];
            if sum(this.param_class == this.PAR_REC_EB) > 0
                for r = 1 : n_rec
                    idx_par = find(this.class_par == this.PAR_REC_EB & this.rec_par == r & ~this.out_par); % one pseudorange bias per reciever
                    if ~isempty(idx_par)
                        
                        % tell which is pseudorange from the electronic bias
                        idx_par_psrange = false(size(idx_par));
                        sys_c_par =  zeros(size(idx_par));
                        for i = 1 : length(idx_par)
                            sys_c_par(i) =  this.unique_obs_codes{this.obs_codes_id_par(idx_par(i))}(1);
                            if this.phase_par(idx_par(i)) == 1
                                idx_par_psrange(i) = true;
                            end
                        end
                        idx_par_phase = idx_par(~idx_par_psrange);
                        
                        % system of the electronic bias
                        sys_c_par_phase = sys_c_par(~idx_par_psrange);
                        sys_c_par_psrange = sys_c_par(idx_par_psrange);
                        idx_par_psrange =  idx_par(idx_par_psrange);
                        
                        
                        if sum(this.param_class == this.PAR_REC_CLK) > 0
                            if ~isempty(idx_par_psrange) % <--- remove one pseudorange because it is less complicated afterwards
                                id_obs_par  = this.obs_codes_id_par(idx_par_psrange);
                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                idx_rm = [idx_rm; uint32(idx_par_psrange(idx_idx_par))]; % <- all bias of the same observation
                                wl_ref = this.wl_id_par(idx_par_psrange(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                sys_c_ref = sys_c_par(idx_idx_par(1));
                                idx_par_psrange(idx_idx_par) = [];
                                sys_c_par_psrange(idx_idx_par) = [];
                            else
                                id_obs_par  = this.obs_codes_id_par(idx_par_phase);
                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                idx_rm = [idx_rm; uint32(idx_par_phase(idx_idx_par))]; % <- all bias of the same observation
                                wl_ref = this.wl_id_par(idx_par_phase(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                sys_c_ref = sys_c_par(idx_idx_par(1));
                                idx_par_phase(idx_idx_par) = [];
                                sys_c_par_phase(idx_idx_par) = [];
                            end
                        end
                        if sum(this.param_class == this.PAR_IONO) > 0 && this.ls_parametrization.iono(2) == LS_Parametrization.SING_REC
                            if ~isempty(idx_par_psrange)
                                for sys_c = unique(sys_c_par_psrange)' % <- each system has its own ionpspherese common to the same elctronic bias so they a new rank deficency is introduced
                                    if sys_c == sys_c_ref
                                        idx_tmp= idx_par_psrange(sys_c_par_psrange == sys_c & wl_ref ~= this.wl_id_par(idx_par_psrange));
                                    else
                                        idx_tmp= idx_par_psrange(sys_c_par_psrange == sys_c);
                                    end
                                    if~isempty(idx_tmp)
                                        id_obs_par  = this.obs_codes_id_par(idx_tmp);
                                        chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                        idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                        idx_rm = [idx_rm; idx_tmp(idx_idx_par)];
                                    end
                                end
                            else
                                for sys_c = unique(sys_c_par_phase)' % <- each system has its own ionpspherese common to the same elctronic bias so they a new rank deficency is introduced
                                    if sys_c == sys_c_ref
                                        idx_tmp= idx_par_phase(sys_c_par_phase == sys_c & wl_ref ~= this.wl_id_par(idx_par_phase));
                                    else
                                        idx_tmp= idx_par_phase(sys_c_par_phase == sys_c);
                                    end
                                    if~isempty(idx_tmp)
                                        id_obs_par  = this.obs_codes_id_par(idx_tmp);
                                        chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                        idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                        idx_rm = [idx_rm; idx_tmp(idx_idx_par)];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            if sum(this.param_class == this.PAR_IONO) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                % for each satellite if there is not ata least one double frequrency receiver observing the satellite remove the iono paramter, the delay is going to be absorbed by the clock
                for s = this.unique_sat_goid
                    idx_par = find(this.class_par == this.PAR_IONO & this.sat_par == s & ~this.out_par); % two code biases not from the same frequency
                    n_freq = numel(unique(this.wl_id_obs(this.satellite_obs == s)));
                    if n_freq == 1
                        idx_rm = [idx_rm; idx_par];
                    end
                end
            end
            
            % remove ones bias per satellite
            if  sum(this.param_class == this.PAR_SAT_EB) > 0
                for s = this.unique_sat_goid
                    idx_par = find(this.class_par == this.PAR_SAT_EB & this.sat_par == s & ~this.out_par); % two code biases not from the same frequency
                    if any(idx_par)
                        idx_par_psrange = false(size(idx_par));
                        for i = 1 : length(idx_par)
                            if this.phase_par(idx_par(i)) == 1
                                idx_par_psrange(i) = true;
                            end
                        end
                        idx_par_phase = idx_par(~idx_par_psrange);
                        idx_par_psrange =  idx_par(idx_par_psrange);
                        if sum(this.param_class == this.PAR_SAT_CLK) > 0
                            if ~isempty(idx_par_psrange) % <--- remove one pseudorange because it is less complicated afterwards
                                id_obs_par  = this.obs_codes_id_par(idx_par_psrange);
                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                idx_rm = [idx_rm; uint32(idx_par_psrange(idx_idx_par))]; % <- all bias of the same observation
                                wl_ref = this.wl_id_par(idx_par_psrange(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                this.log.addMessage(this.log.indent(sprintf('Pseudorange %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_psrange(idx_idx_par(1)))},s)));
                                %idx_par_psrange(idx_idx_par) = [];
                            else
                                id_obs_par  = this.obs_codes_id_par(idx_par_phase);
                                chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                idx_idx_par = find(id_obs_par == chosen_id_obs);
                                idx_rm = [idx_rm; uint32(idx_par_phase(idx_idx_par))]; % <- all bias of the same observation
                                wl_ref = this.wl_id_par(idx_par_phase(idx_idx_par(1))); % <- you can not then remove from  the same frequency and system
                                sys_c_ref = sys_c_par(idx_idx_par(1));
                                %idx_par_phase(idx_idx_par) = [];
                                this.log.addMessage(this.log.indent(sprintf('Phase %s choosen as reference for sat %d',wl_ref,s)));
                            end
                        end
                        if sum(this.param_class == this.PAR_IONO) > 0
                            if ~isempty(idx_par_psrange)
                                idx_par_psrange = idx_par_psrange(wl_ref ~= this.wl_id_par(idx_par_psrange)); %<- remove bias of the same frequency of one thta has been already removed
                                if ~isempty(idx_par_psrange)
                                    id_obs_par  = this.obs_codes_id_par(idx_par_psrange);
                                    chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                    idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                    idx_rm = [idx_rm; uint32(idx_par_psrange(idx_idx_par))];
                                    this.log.addMessage(this.log.indent(sprintf('Pseudorange %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_psrange(idx_idx_par(1)))},s)));
                                end
                            else
                                idx_par_phase = idx_par_phase(wl_ref ~= this.wl_id_par(idx_par_phase));  %<- remove bias of the same frequency of one thta has been already removed
                                if ~isempty(idx_par_phase)
                                    id_obs_par  = this.obs_codes_id_par(idx_par_phase);
                                    chosen_id_obs = mode(id_obs_par(id_obs_par~=0));
                                    idx_idx_par = find(id_obs_par == chosen_id_obs);  % <- all bias of the same observation
                                    idx_rm = [idx_rm; uint32(idx_par_phase(idx_idx_par))];
                                    this.log.addMessage(this.log.indent(sprintf('Phase %s choosen as reference for sat %d',this.unique_obs_codes{this.obs_codes_id_par(idx_par_phase(idx_idx_par(1)))},s)));
                                end
                            end
                        end
                    end
                end
            end
            
            if  sum(this.param_class == this.PAR_SAT_EB) > 0 && sum(this.param_class == this.PAR_REC_EB) > 0
                idx_par = find(this.class_par == this.PAR_REC_EB & ~this.out_par);
                rec_wl_id = uint64(this.wl_id_par(idx_par)) + 100*uint64(this.rec_par(idx_par));
                rec_wl_id_rm = uint64(this.wl_id_par(idx_rm)) + 100*uint64(this.rec_par(idx_rm));
                [conflictig_id] = ismember(rec_wl_id,rec_wl_id_rm); % <- you can not remove from same receiver and frequency this could generate inconsistency
                idx_par(conflictig_id) = [];
                if ~isempty(idx_par)
                    idx_rm = [idx_rm; uint32(idx_par(1))];
                    this.log.addMessage(this.log.indent(sprintf('Receiver %d and %s choosen as reference',this.rec_par(idx_par(1)),this.unique_obs_codes{this.obs_codes_id_par(idx_par(1))})));
                end
            end
            
            
            
            
            
            
            % remove one bias per signal from one recievr
            if  sum(this.param_class == this.PAR_SAT_EB) > 0 && sum(this.param_class == this.PAR_REC_EB) > 0
                for e = 1: length(this.unique_obs_codes)
                    if sum(this.obs_codes_id_par(idx_rm(this.class_par(idx_rm) == this.PAR_REC_EB)) == e) == 0% <- if it has not been previously removed
                        idx_par = find(this.class_par == this.PAR_REC_EB & this.obs_codes_id_par == e & ~this.out_par);
                        if ~isempty(idx_par)
                            idx_rm = [idx_rm; uint32(idx_par(1))];
                            this.log.addMessage(this.log.indent(sprintf('Receiver %d choosen as reference for obs %s',this.rec_par(idx_par(1)),this.unique_obs_codes{e})));
                        end
                    end
                end
            end
            
            % remove one linear trend bias per signal from one recievr
            if  sum(this.param_class == this.PAR_SAT_EB) > 0 && sum(this.param_class == this.PAR_REC_EB_LIN) > 0
                for e = 1: length(this.unique_obs_codes)
                    idx_par = find(this.class_par == this.PAR_REC_EB_LIN & this.obs_codes_id_par == e & ~this.out_par);
                    if ~isempty(idx_par)
                        idx_rm = [idx_rm; uint32(idx_par(1))];
                        this.log.addMessage(this.log.indent(sprintf('Receiver %d choosen as reference for obs %s',this.rec_par(idx_par(1)),this.unique_obs_codes{e})));
                    end
                end
            end
            
            % remove one clock per epoch
            u_ep = unique(this.time_par);
            if sum(this.param_class == this.PAR_REC_CLK) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                for e = u_ep'
                    idx_par = find(this.class_par == this.PAR_REC_CLK & this.time_par(:,1) == e & ~this.out_par);
                    if length(idx_par)>1
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                    end
                end
            end
            % phase only TBD!!
            
            
            % for each sat and fro each contiguous set of ambiguity remove
            % one abiguity per set of phase bias
            
            % remember -> there is the possibility of island of
            % disconnected tracking that you have not considered
            
            %rule -> you can not have two full receiver removed at this
            %stage
            % you start from the most populus receiver than you move to the
            % other keeping in mind that after the first no receiver can be
            % completely removed
            if sum(this.param_class == this.PAR_AMB) > 0 && (sum(this.param_class == this.PAR_SAT_CLK) > 0 || sum(this.param_class == this.PAR_REC_CLK) > 0)
                idx_amb_rm_sat = [];
                if (sum(this.param_class == this.PAR_SAT_CLK) > 0)
                    % find the elecronic bias assoictaed with each ambiguity
                    idx_ambs = find(this.class_par == this.PAR_AMB);
                    amb2eb = zeros(size(idx_ambs));
                    % find to which electronic bias the ambiguity is tied
                    for e = 1: length(idx_ambs)
                        idx_obs_sample = find(this.A_idx(:,this.param_class == this.PAR_AMB) == idx_ambs(e),1,'first');
                        ebs_tmp = this.obs_codes_id_par(this.A_idx(idx_obs_sample, this.param_class == this.PAR_SAT_EB));
                        amb2eb(e) = ebs_tmp(1);
                    end
                    clearvars ebs_tmp
                    sat_eb_const =  this.ls_parametrization.sat_eb(1) == LS_Parametrization.CONST;
                    if sat_eb_const
                        jmps_sat ={};
                        jmps_sat_el ={}; % have been elimated an ambiguity from the block
                        % determine all arcs jum
                        for s = this.unique_sat_goid
                            idx_par = this.class_par == this.PAR_AMB & this.sat_par ==  s & ~this.out_par;
                            idx_par = find(idx_par);
                            time_par = this.time_par(idx_par,:);
                            sat_par  = this.sat_par(idx_par,:);
                            obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                            u_time = unique(time_par);
                            amb2arc_a = 1000*uint32(sat_par) + uint32(obs_codes_id_par);
                            u_arc = unique(amb2arc_a);
                            arc2eb = rem(u_arc,1000);
                            
                            
                            [~,amb2arc] = ismember(amb2arc_a,u_arc);
                            u_eb = unique(arc2eb);
                            eb_arc_rem = false(size(u_eb));
                            amb_mat = zeros(max(u_time),length(u_arc),'uint32');
                            for t = 1 : size(time_par,1)
                                amb_mat(time_par(t,1)+1:time_par(t,2),amb2arc(t)) = idx_par(t);
                            end
                            jmps = [(find(diff(sum(amb_mat,2) > 0) == 1) +1); max(u_time)];
                            if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                                jmps = [1; jmps];
                            end
                            jmps_sat{s} = jmps;
                            jmps_sat_el{s} = false(length(jmps)-1,1);
                            
                        end
                    end
                    ebs = unique(amb2eb)';
                    for eb = ebs
                        if sum(this.param_class == this.PAR_AMB) > 0 && sum(this.param_class == this.PAR_SAT_CLK) > 0
                            is_first_complete = false; %flag to know if the recievr has been removed completely
                            idx_ambs_e = idx_ambs;
                            idx_ambs_e(amb2eb ~= eb) = [];
                            rec_sat_mtx = zeros(n_rec,n_sat);
                            for aa = idx_ambs_e
                                rec_sat_mtx(this.rec_par(aa),this.sat_par(aa)) = 1;
                            end
                            idx_ambs_e = Core_Utils.ordinal2logical(idx_ambs_e,length(this.class_par));
                            
                            [~, rec_preference] = sort(sum(rec_sat_mtx,2),'descend');
                            for s = this.unique_sat_goid
                                
                                idx_par = idx_ambs_e & this.sat_par == s & ~this.out_par;
                                idx_par = find(idx_par);
                                if any(idx_par)
                                    time_par = this.time_par(idx_par,:);
                                    rec_par  = this.rec_par(idx_par,:);
                                    obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                                    
                                    u_time = unique(time_par);
                                    u_rec = unique(rec_par);
                                    amb2arc = 1000*uint32(rec_par) + uint32(obs_codes_id_par);
                                    u_arc = unique(amb2arc);
                                    arc2eb = rem(u_arc,1000);
                                    
                                    
                                    [~,amb2u_arc] = ismember(amb2arc,u_arc);
                                    u_eb = unique(arc2eb);
                                    
                                    amb_mat = zeros(max(u_time),length(u_arc),'uint32');
                                    rec_amb_mat =  zeros(1,length(u_arc),'uint32');
                                    for t = 1 : size(time_par,1)
                                        amb_mat(time_par(t,1)+1:time_par(t,2),amb2u_arc(t)) = idx_par(t);
                                        rec_amb_mat(amb2u_arc(t)) = floor(amb2arc(t) / 1000);
                                    end
                                    if ~sat_eb_const
                                        jmps = [(find(diff(sum(amb_mat,2) > 0) == 1) +1); max(u_time)];
                                        if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                                            jmps = [1; jmps];
                                        end
                                    else
                                        if any(~jmps_sat_el{s})
                                            jmps = jmps_sat{s};
                                        else
                                            jmps = [1; size(amb_mat,1)];
                                        end
                                    end
                                    
                                    
                                    for j = 1 : (length(jmps) -1)
                                        first_rem = true; % one has to be removed al the times
                                        
                                        jmp_s = jmps(j);
                                        jmp_e = jmps(j+1);
                                        ambs = amb_mat(jmp_s:min(size(amb_mat,1),jmp_e),:);
                                        if any(any(ambs))
                                            rr = 1;
                                            not_found = true;
                                            while rr <= length(rec_preference) && not_found
                                                rp = rec_preference(rr);
                                                if ~(sum(rec_sat_mtx(rec_preference(1),:) == 1) == 0 && sum(rec_sat_mtx(rp,:) == 1) <= 1)
                                                    ambs_r = ambs(:,rec_amb_mat == rp);
                                                    if any(any(ambs_r)) && (~sat_eb_const || ~jmps_sat_el{s}(j) ||  first_rem )
                                                        idx_poss_amb = mode(noZero(ambs_r(:)));
                                                        idx_amb_rm_sat = [idx_amb_rm_sat; uint32(idx_poss_amb)];
                                                        if sat_eb_const
                                                            jmps_sat_el{s}(j) = true;
                                                        end
                                                        first_rem = false;
                                                        if ~sat_eb_const || ~any(~jmps_sat_el{s})
                                                            rec_sat_mtx(rr,s) == 2; % two means eliminated
                                                        end
                                                        not_found = false;
                                                    end
                                                end
                                                rr = rr +1;
                                            end
                                        end
                                    end
                                    
                                end
                                idx_rm = [idx_rm; uint32(idx_amb_rm_sat)];
                            end
                        end
                    end
                end
                if (sum(this.param_class == this.PAR_REC_CLK) > 0)
                    % find the elecronic bias assoictaed with each ambiguity
                    idx_ambs = find(this.class_par == this.PAR_AMB);
                    amb2eb = zeros(size(idx_ambs));
                    % find to which electrinuc bias the ambiguity is tied
                    for e = 1: length(idx_ambs)
                        idx_obs_sample = find(this.A_idx(:,this.param_class == this.PAR_AMB) == idx_ambs(e),1,'first');
                        amb2eb(e) = this.obs_codes_id_par(this.A_idx(idx_obs_sample,this.param_class == this.PAR_REC_EB));
                    end
                    ebs = unique(amb2eb)';
                    rec_eb_const =  this.ls_parametrization.rec_eb(1) == LS_Parametrization.CONST;
                    if rec_eb_const
                        jmps_rec ={};
                        jmps_rec_el ={}; % have been elimated an ambiguity from the block
                        % determine all arcs jum
                        for r = 1: size(this.rec_xyz,1);
                            idx_par = this.class_par == this.PAR_AMB & this.rec_par ==  r & ~this.out_par;
                            idx_par = find(idx_par);
                            time_par = this.time_par(idx_par,:);
                            sat_par  = this.sat_par(idx_par,:);
                            obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                            u_time = unique(time_par);
                            amb2arc_a = 1000*uint32(sat_par) + uint32(obs_codes_id_par);
                            u_arc = unique(amb2arc_a);
                            arc2eb = rem(u_arc,1000);
                            
                            
                            [~,amb2arc] = ismember(amb2arc_a,u_arc);
                            u_eb = unique(arc2eb);
                            eb_arc_rem = false(size(u_eb));
                            amb_mat = zeros(max(u_time),length(u_arc),'uint32');
                            for t = 1 : size(time_par,1)
                                amb_mat(time_par(t,1)+1:time_par(t,2),amb2arc(t)) = idx_par(t);
                            end
                            jmps = [(find(diff(sum(amb_mat,2) > 0) == 1) +1); max(u_time)];
                            if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                                jmps = [1; jmps];
                            end
                            jmps_rec{r} = jmps;
                            jmps_rec_el{r} = false(length(jmps)-1,1);
                            
                        end
                    end
                    for eb = ebs
                        % for each rec and for each contiguos set of ambiguity remove one
                        if sum(this.param_class == this.PAR_AMB) > 0 && sum(this.param_class == this.PAR_REC_CLK) > 0
                            o_ch = this.obs_codes_id_par(idx_amb_rm_sat);
                            o_ch(o_ch < 0) = length(this.unique_obs_codes) - o_ch(o_ch < 0); % negative index stand for set  put them positive after the numebr of type of observation
                            forbidden_arc = 1e6*uint32(this.rec_par(idx_amb_rm_sat)) + 1000*uint32(this.sat_par(idx_amb_rm_sat)) + uint32(this.ls_parametrization.rec_eb(4) == LS_Parametrization.SING_TRACK)*uint32(o_ch); % this ambiguities can not be elemitaing without genrating "tension" in the system  %+ uint32(this.obs_codes_id_par(idx_amb_rm_sat)) if there is a bias per tracking the condition is weaker
                            for r = 1: size(this.rec_xyz,1);
                                forbidden_arc_rec = rem(forbidden_arc(floor(forbidden_arc/1e6) == r),1e6);
                                idx_ambs_e = idx_ambs;
                                idx_ambs_e(amb2eb ~= eb) = [];
                                idx_ambs_e = Core_Utils.ordinal2logical(idx_ambs_e,length(this.class_par));
                                idx_par = idx_ambs_e & this.rec_par ==  r & ~this.out_par;
                                idx_par = find(idx_par);
                                if any(idx_par)
                                    
                                    time_par = this.time_par(idx_par,:);
                                    sat_par  = this.sat_par(idx_par,:);
                                    obs_codes_id_par = this.obs_codes_id_par(idx_par,:);
                                    u_time = unique(time_par);
                                    amb2arc_a = 1000*uint32(sat_par) + uint32(obs_codes_id_par);
                                    u_arc = unique(amb2arc_a);
                                    arc2eb = rem(u_arc,1000);
                                    
                                    
                                    [~,amb2arc] = ismember(amb2arc_a,u_arc);
                                    amb_mat = zeros(max(u_time),length(u_arc),'uint32');
                                    for t = 1 : size(time_par,1)
                                        amb_mat(time_par(t,1)+1:time_par(t,2),amb2arc(t)) = idx_par(t);
                                    end
                                    if ~rec_eb_const
                                        jmps = [(find(diff(sum(amb_mat,2) > 0) == 1) +1); max(u_time)];
                                        if ~isempty(amb_mat) && sum(abs(amb_mat(1,:) )) ~= 0 %<- if first epoch is full start of the arc is not detected
                                            jmps = [1; jmps];
                                        end
                                    else
                                        if any(~jmps_rec_el{r})
                                            jmps = jmps_rec{r};
                                        else
                                            jmps = [1; size(amb_mat,1)];
                                        end
                                    end
                                    for j = 1 : (length(jmps) -1)
                                        first_rem = true; % one has to be removed al the times
                                        jmp_s = jmps(j);
                                        jmp_e = jmps(j+1);
                                        ambs = amb_mat(jmp_s:min(size(amb_mat,1),jmp_e),:);
                                        if any(any(ambs))
                                            [id_poss_rm]  = mode(noZero(ambs(:)));
                                            idx_start = sum(ambs == id_poss_rm) > 0; % i  case everything has been removed to exit the loop
                                            while any(any(ambs)) && sum(id_poss_rm ==  idx_rm) > 0 | sum(amb2arc_a(find(idx_par == id_poss_rm)) == forbidden_arc_rec) > 0 ...
                                                    | (this.ls_parametrization.rec_eb(4) ~= LS_Parametrization.SING_TRACK && sum(floor(amb2arc_a(find(idx_par == id_poss_rm))/1000) == floor(forbidden_arc_rec/1000)) > 0)% it might be that the ambiguity was previouly removed in the satellite round
                                                ambs(ambs == id_poss_rm) = 0;
                                                id_poss_rm = mode(noZero(ambs(:)));
                                            end
                                            if id_poss_rm > 0 && (~rec_eb_const || ~jmps_rec_el{r}(j) ||  first_rem )
                                                idx_rm = [idx_rm; uint32(id_poss_rm)];
                                                if rec_eb_const
                                                    jmps_rec_el{r}(j) = true;
                                                end
                                                first_rem = false;
                                                
                                            end
                                            
                                            
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % ---- (multi receiver) for each epoche remove one coordinate --------
            if (sum(this.param_class == this.PAR_REC_X) > 0 || sum(this.param_class == this.PAR_REC_Y) > 0  ||  sum(this.param_class == this.PAR_REC_Z) > 0 || sum(this.param_class == this.PAR_TROPO) > 0  || sum(this.param_class == this.PAR_TROPO_E) > 0  || sum(this.param_class == this.PAR_TROPO_N) > 0) && sum(this.param_class == this.PAR_SAT_CLK) > 0
                for e = u_ep'
                    idx_e = this.time_par(:,1) >= e & ( this.time_par(:,1)==0 | this.time_par(:,2) < e) ;
                    idx_par = this.class_par == this.PAR_REC_X & idx_e & ~this.out_par;
                    idx_par = find(idx_par);
                    if any(any(idx_par))
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                    idx_par = this.class_par == this.PAR_REC_Y & idx_e & ~this.out_par;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; idx_par(idx_rm_rm)];
                        end
                    end
                    idx_par = this.class_par == this.PAR_REC_Z & idx_e & ~this.out_par;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                        end
                    end
                    idx_par = this.class_par == this.PAR_TROPO & idx_e & ~this.out_par;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                        end
                    end
                    idx_par = this.class_par == this.PAR_TROPO_E & idx_e & ~this.out_par;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                        end
                    end
                    idx_par = this.class_par == this.PAR_TROPO_N & idx_e & ~this.out_par;
                    idx_par = find(idx_par);
                    
                    if any(idx_par)
                        [~,idx_rm_rm] = min(this.rec_par(idx_par));
                        if sum(idx_par(idx_rm_rm) == idx_rm) == 0
                            idx_rm = [idx_rm; uint32(idx_par(idx_rm_rm))];
                        end
                    end
                end
            end
            this.idx_rd = unique(noZero(idx_rm));
            %             for i = length(this.idx_rd) : -1 : 1
            %                 ir = this.idx_rd(i);
            %                 this.A_idx(this.A_idx == ir) = 0;
            %                 i_mj = this.A_idx > ir;
            %                 this.A_idx(i_mj) = this.A_idx(i_mj) -1;
            %             end
        end
        
        
        function absValRegularization(this,p_class, var)
            % regularize parameters to zero (Tykhnov aka ridge aka L2)
            %
            % this.absValRegularization(param_id, var)
            par_ids = this.param_class == p_class;
            u_p_id = unique(this.A_idx(:,par_ids));
            n_par = length(u_p_id);
            A_idx_tmp = zeros(n_par, 2,'uint32'); % tykhonv regualrization are now limited to two parameters
            A_tmp = zeros(n_par, 2);
            A_tmp(:,1) = 1;
            A_idx_tmp(:,1) = u_p_id;
            this.A_pseudo = [this.A_pseudo; A_tmp];
            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
            this.variance_pseudo = [this.variance_pseudo; var*ones(n_par,1)];
            this.receiver_pseudo = [this.receiver_pseudo; this.rec_par(u_p_id)];
            if isempty(this.time_pseudo)
                this.time_pseudo = GPS_Time(this.time_min.getMatlabTime + this.time_par(u_p_id));
            else
                this.time_pseudo.addEpoch(GPS_Time(this.time_min.getMatlabTime + this.time_par(u_p_id)));
            end
            this.satellite_pseudo = [this.satellite_pseudo; this.sat_par(u_p_id)];
        end
        
        function timeRegularization(this, p_class, var_per_sec)
            % first order tykhonv regualrization in time
            %
            % this.timeRegularization(this, param_id, var)
            par_ids = this.param_class == p_class;
            p_idx = unique(find(this.class_par == p_class));
            rec_idx = this.rec_par(p_idx);
            u_rec = unique(rec_idx);
            sat_idx = this.sat_par(p_idx);
            u_sat = unique(sat_idx);
            ch_idx  =  this.obs_codes_id_par(p_idx);
            u_ch = unique(ch_idx);
            for r = u_rec'
                for s = u_sat'
                    for c = u_ch'
                        p_tmp = p_idx(rec_idx == r & sat_idx == s & ch_idx == c); % find the idx of the obsrevations
                        if ~isempty(p_tmp)
                            u_p_id = unique(p_tmp);
                            n_par = length(u_p_id);
                            A_tmp = [ones(n_par-1,1) -ones(n_par-1,1)];
                            A_idx_tmp = [u_p_id(1:(end-1)) u_p_id(2:end)];
                            this.A_pseudo = [this.A_pseudo; A_tmp];
                            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
                            
                            this.variance_pseudo = [this.variance_pseudo; var_per_sec * diff(double(this.time_par(u_p_id)))];
                            % taking the indices of first epoch
                            this.receiver_pseudo = [this.receiver_pseudo; this.rec_par(u_p_id(1:end-1))];
                            if isempty( this.time_pseudo)
                                this.time_pseudo = GPS_Time(this.time_min.getMatlabTime + this.time_par(u_p_id(1:end-1)));
                            else
                                this.time_pseudo.addEpoch(this.time_min.getMatlabTime + this.time_par(u_p_id(1:end-1)));
                            end
                            this.satellite_pseudo = [this.satellite_pseudo; this.sat_par(u_p_id(1:end-1))];
                        end
                    end
                end
            end
        end
        
        function spatialRegularization(this, law)
            % tykhonv regualrization in space
            %
            % this.spatialRegularization(this, law)
            
            
            n_rec_tot = size(this.xyz,1);
            dist_mat = zeros(n_rec_tot);
            for i = 1 : n_rec_tot
                for j = (i+1) : n_rec_tot
                    coord1 = Coordinates();
                    coord1.xyz = this.xyz(i,:);
                    coord2 = Coordinates();
                    coord2.xyz = this.xyz(j,:);
                    dist_mat(i,j) = coord1.ellDistanceTo(coord2);
                end
            end
            var_mat = law(dis_mat);
            
            % -- construct A and indices
            par_ids = this.parameter_class == p_class;
            p_idx = unique(this.A_idx(:,par_ids));
            sat_idx = this.sat_par(p_idx);
            u_sat = unique(sat_idx);
            ch_idx  =  this.ch_par(p_idx);
            time_idx = this.time_par(p_idx);
            u_ch = unique(ch_idx);
            for s = u_sat
                for c = u_ch
                    p_tmp = p_idx(rec_idx == r & sta_idx == s & ch_idx == c); % find the idx of the obsrevations
                    u_ep = unique(this.time_par(p_tmp));
                    for e = u_ep
                        p_tmp = p_idx(rec_idx == r & sta_idx == s & ch_idx == c & time_idx == e); % find the idx of the obsrevations
                        if ~isempty(p_tmp)
                            u_p_id = unique(p_tmp);
                            n_par = length(u_p_id);
                            rec_par = this.rec_par(u_p_id);
                            n_rec = length(rec_apr);
                            n_pseudobs = n_rec*(n_rec - 1);
                            A_tmp = [ones(n_par-1,1) -ones(n_par-1,1)];
                            idx_1 = [];
                            for r = 1 : n_rec
                                idx_1 = [idx_1; (r:n_rec)'];
                            end
                            idx_2 = idx_1 + 1;
                            A_idx_tmp = [u_p_id(idx_1) u_p_id(idx_2)];
                            var = var_mat(rec_par(idx_1)+n_rec_tot*(rec_par(idx_1) -1));
                            this.A_pseudo = [this.A_pseudo; A_tmp];
                            this.A_idx_pseudo = [this.A_idx_pseudo; A_idx_tmp];
                            
                            this.variance_pseudo = [this.variance_pseudo; var_per_sec * var];
                            % taking the indices of first epoch
                            this.receiver_pseudo = [this.receiver_pseudo; zeros(n_pseudo_obs,1)];
                            this.time_pseudo.addEpoch(GPS_Time(this.time_min + e*ones(n_pseudo_obs,1)));
                            this.satellite_pseudo = [this.satellite_pseudo; s*ones(n_pseudo_obs,1)];
                        end
                    end
                end
            end
        end
        
        function hemisphereRegularization(this, law)
            % tykhonv regualrization on the hemisphere
            %
            % this.spatialRegularization(this, law)
            % TBD
        end
        
        function reduceForNuisanceParameters(this, param_id)
            % reduce for the parameter
        end
        
        
        function solve(this, fix)
            % %solve the least squares
            % SYNTAX:
            %    this.solve(<fix>)
            % INPUT:
            %    fix : seek to fix integer parameters
            if nargin < 2
                fix = false;
            end
            % ------ mark short arc as outlier
            this.markShortArcs(1);
            % ------ mark single reciver or satellite epoch as outlier
            this.markSingledObs();
            % ------ remove outlier
            this.generateOutliedParIdx();
            % ------ remove full rank deficencies
            this.removeFullRankDeficency();
            % ------ form the normal matrix
            n_obs = size(this.A,1) + size(this.A_pseudo,1);
            n_par = double(max(max(this.A_idx)));
            rows = repmat((1:size(this.A,1))',1,size(this.A,2));
            rows_pseudo = repmat(size(this.A,1)+(1:size(this.A_pseudo,1))',1,size(this.A_pseudo,2));
            rows = [rows(:); rows_pseudo(:)];
            columns = double(zero2n([this.A_idx(:); this.A_idx_pseudo(:)],1));
            values = [this.A(:); this.A_pseudo(:)];
            A = sparse(rows, columns, values, n_obs, n_par);
            %this.idx_rd = [];
            n_out = sum(this.outlier_obs);
            A_out = A(this.outlier_obs > 0,:);
            A(:, [this.idx_rd; find(this.out_par)]) = [];
            
            A(this.outlier_obs > 0, :) = [];
            
            class_par = this.class_par;
            class_par([this.idx_rd; find(this.out_par)]) = [];
            
            zero_pars = sum(A~=0) == 0;
            A(:,zero_pars) = [];
            class_par(zero_pars) = [];
            valid_pars = find(~Core_Utils.ordinal2logical([this.idx_rd; find(this.out_par)],n_par)); % sometimes with splines spme paramter have a zero entry
            this.idx_rd = [this.idx_rd; valid_pars(zero_pars)];
            vars = [1./this.variance_obs(~this.outlier_obs); 1./this.variance_pseudo];
            mean_vars = 1 ; %mean(vars);
            vars = vars ./ mean_vars;
            Cyy =  spdiags(vars,0,n_obs - n_out,n_obs - n_out);
            x_est = zeros(n_par -length(this.idx_rd) - sum(this.out_par),1);
            Aw = A'*Cyy;
            N = Aw*A;
            y = sparse([this.obs(~this.outlier_obs); zeros(size(this.A_pseudo,1),1)]);
            B = Aw*y;
            
            clearvars Aw
            % ------ reduce for sat clock, rec clock and iono
            idx_reduce_sat_clk = class_par == this.PAR_SAT_CLK;
            idx_reduce_rec_clk = class_par == this.PAR_REC_CLK;
            idx_reduce_iono = class_par == this.PAR_IONO;
            
            iono = sum(class_par  == this.PAR_IONO) > 0;
            if iono
                n_iono = sum(idx_reduce_iono);
                iIono = spdiags(1./diag(N(idx_reduce_iono,idx_reduce_iono)),0,n_iono,n_iono);
                Nx_iono = N(~idx_reduce_iono,idx_reduce_iono); % cross term reduce iono
                Nt = Nx_iono * iIono;
                N = N(~idx_reduce_iono,~idx_reduce_iono) - Nt * N(idx_reduce_iono,~idx_reduce_iono);
                B_iono =  B(idx_reduce_iono);
                B = B(~idx_reduce_iono) - Nt * B_iono;
            end
            
            sat_clk = sum(class_par  == this.PAR_SAT_CLK) > 0;
            if sat_clk
                i_sat_clk_tmp = idx_reduce_sat_clk(~idx_reduce_iono);
                n_clk_sat = sum(i_sat_clk_tmp);
                iSatClk = spdiags(1./diag(N(i_sat_clk_tmp,i_sat_clk_tmp)),0,n_clk_sat,n_clk_sat);
                Nx_satclk = N(~i_sat_clk_tmp, i_sat_clk_tmp);
                Nt = Nx_satclk * iSatClk;
                N = N(~i_sat_clk_tmp,~i_sat_clk_tmp) - Nt * N(i_sat_clk_tmp, ~i_sat_clk_tmp);
                B_satclk =  B(i_sat_clk_tmp);
                B = B(~i_sat_clk_tmp) - Nt * B_satclk;
            end
            
            rec_clk = sum(class_par  == this.PAR_REC_CLK) > 0;
            if rec_clk
                i_rec_clk_tmp = idx_reduce_rec_clk(~idx_reduce_iono & ~idx_reduce_sat_clk);
                iRecClk = inv(N(i_rec_clk_tmp,i_rec_clk_tmp));
                Nx_recclk = N(~i_rec_clk_tmp, i_rec_clk_tmp);
                Nt = Nx_recclk * iRecClk;
                N = N(~i_rec_clk_tmp, ~i_rec_clk_tmp) - Nt * N(i_rec_clk_tmp, ~i_rec_clk_tmp);
                B_recclk = B(i_rec_clk_tmp);
                B = B(~i_rec_clk_tmp) - Nt * B_recclk;
            end
            
            x_reduced = N\B;
            
            % ------- fix the ambiguities
            
            if sum(this.param_class == this.PAR_AMB) > 0 && fix || true
                % get the ambiguity inverse matrxi
                idx_amb = find(class_par(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono) == this.PAR_AMB);
                if any(idx_amb)
                    n_amb = length(idx_amb);
                    amb_y  = sparse(idx_amb,1:length(idx_amb),ones(size(idx_amb)),size(N,1),numel(idx_amb));
                    C_amb_amb = N \ amb_y;
                    idx_rm_line = true(size(C_amb_amb,1),1);
                    idx_rm_line(idx_amb) = false;
                    C_amb_amb(idx_rm_line,:) = [];
                    amb_float = x_reduced(idx_amb);
                    amb_fixed = amb_float;
                    l_fixed = abs(fracFNI(amb_float)) < 0.1 & diag(C_amb_amb) < 0.01;
                    amb_fixed(l_fixed) = round(amb_fixed(l_fixed));
                    is_fixed = true;
                    %[amb_fixed, is_fixed, l_fixed] = Fixer.fix(full(amb_float), full(C_amb_amb), 'lambda_ILS' );
                    flag_debug = false;
                    if flag_debug
                        N_old = N;
                        B_old = B;
                        figure;
                        coo =x_reduced(1:3);
                        
                        plot3(coo(1),coo(2),coo(3),'o')
                        grid on;
                        hold on;
                        
                        
                        
                        for a = 1 : size(amb_fixed,2)
                            idx_amb_fixed = idx_amb(l_fixed(:, a));
                            n_amb_fixed = sum(l_fixed(:, a));
                            idx_not_amb_fixed = true(size(N,1), 1);
                            idx_not_amb_fixed(idx_amb_fixed) = false;
                            Nt = N(idx_not_amb_fixed, idx_amb_fixed);
                            B(idx_amb_fixed) = [];
                            B = B - sum(Nt * spdiags(amb_fixed(l_fixed(:,a),a),0,n_amb_fixed,n_amb_fixed),2);
                            N(idx_amb_fixed, :) = [];
                            N(:, idx_amb_fixed) = [];
                            
                            x_fixed = N \ B;
                            coo =x_fixed(1:3);
                            plot3(coo(1),coo(2),coo(3),'.','MarkerSize',10)
                            text(coo(1)+0.01,coo(2)+0.01,coo(3)+0.01,num2str(a));
                            N = N_old;
                            B = B_old;
                        end
                        
                    end
                    if is_fixed
                        idx_amb_fixed = idx_amb(l_fixed(:, 1));
                        n_amb_fixed = sum(l_fixed(:, 1));
                        idx_not_amb_fixed = true(size(N,1), 1);
                        idx_not_amb_fixed(idx_amb_fixed) = false;
                        Nt = N(idx_not_amb_fixed, idx_amb_fixed);
                        B(idx_amb_fixed) = [];
                        B = B - sum(Nt * spdiags(amb_fixed(l_fixed(:,1),1),0,n_amb_fixed,n_amb_fixed),2);
                        N(idx_amb_fixed, :) = [];
                        N(:, idx_amb_fixed) = [];
                        
                        x_fixed = N \ B;
                        
                        x_reduced(idx_amb_fixed) = amb_fixed(l_fixed(:,1),1);
                        x_reduced(idx_not_amb_fixed) = x_fixed;
                        clearvars x_fixed Nt C_amb_amb N B
                    end
                end
                
            end
            
            % ------- substitute back
            x_est(~idx_reduce_sat_clk & ~idx_reduce_rec_clk & ~idx_reduce_iono) = x_reduced;
            
            % receiver clock
            if rec_clk
                B_recclk = B_recclk - sum(Nx_recclk' * spdiags(x_reduced,0,length(x_reduced),length(x_reduced)),2);
                x_rec_clk = iRecClk * B_recclk;
                x_est(idx_reduce_rec_clk) = x_rec_clk;
            end
            
            % satellite clcok
            if sat_clk
                n_sat_clk = size(B_recclk,1);
                idx_est = ~ idx_reduce_iono & ~idx_reduce_sat_clk ;
                B_satclk = B_satclk -   sum(Nx_satclk' * spdiags(x_est(idx_est),0,sum(idx_est),sum(idx_est)),2);
                x_sat_clk = iSatClk * B_satclk;
                x_est(idx_reduce_sat_clk) = x_sat_clk;
            end
            
            % iono
            if iono
                n_iono = size(B_iono,1);
                idx_est = ~ idx_reduce_iono;
                B_iono = B_iono -   sum(Nx_iono' * spdiags(x_est(idx_est),0,sum(idx_est),sum(idx_est)),2);
                x_iono = iIono * B_iono;
                x_est(idx_reduce_iono) = x_iono;
            end
            x = zeros(n_par,1);
            idx_est = true(n_par,1);
            idx_est([this.idx_rd ; find(this.out_par)]) = false;
            x(idx_est) = x_est;
            res = nan(size(this.obs));
            
            % generate esatimations also for the out par (to get a residual)
            if n_out > 0
                res_out = this.obs(this.outlier_obs) - A_out(:,~this.out_par & ~Core_Utils.ordinal2logical(this.idx_rd,n_par))*x_est;
                red_out = res_out;
                A_res_red = A_out(:,this.out_par);
                idx_empty = sum(A_res_red,2) == 0;
                red_out(idx_empty) = [];
                A_res_red(idx_empty,:) = [];
                x_out = A_res_red \ red_out;
                res_out(~idx_empty) = red_out - A_res_red * x_out;
                x(this.out_par) = x_out;
                res(this.outlier_obs) = res_out;
            end
            
            res(~this.outlier_obs) = this.obs(~this.outlier_obs) - A(1:sum(~this.outlier_obs),:)*x_est;
            this.res = res;
            this.x = x;
            
        end
        
        function applyWeightingStrategy()
        end
        
        function res = getResidual(this)
        end
        
        function simpleSnoop(this, ph_thr, pr_thr)
            % simple threshold on residual
            %
            % this.Snoop(this, ph_thr, pr_thr)
            idx_out_ph = this.phase_obs & abs(this.res) > ph_thr;
            this.outlier_obs(idx_out_ph) = true;
            idx_out_pr = this.phase_obs == 0 & abs(this.res) > pr_thr;
            this.outlier_obs(idx_out_pr) = true;
        end
        
        function [res_ph, sat, obs_id, res_id] = getPhRes(this, rec)
            % Get phase residual
            %
            % OUPUT
            %   res_ph          matrix of phase residuals
            %   sat             go_id of the satellite
            %   obs_id          id of the array ls.unique_obs_codes indicating the constallation and tracking of the column
            %   id_res          indix of the value in the res array
            %
            % SYNTAX
            %   [res_ph, sat, obs_id] = this.getPhRes(rec)
            if nargin <2
                rec = 1;
            end
            idx_rec = find(this.receiver_obs == rec);
            u_stream = unique(uint32(this.satellite_obs(idx_rec  & this.phase_obs )) + 1000*uint32(this.obs_codes_id_obs(idx_rec  & this.phase_obs )));
            n_stream = length(u_stream);
            time_res = this.time_obs.getNominalTime.getEpoch(idx_rec).minimum;
            duration = this.time_obs.getNominalTime.getEpoch(idx_rec).maximum - time_res;
            time_res.addSeconds(0:this.time_obs.getRate:duration);
            res_ph = nan(time_res.length,n_stream);
            res_id = zeros(time_res.length,n_stream,'uint32');
            sat = nan(1,n_stream);
            obs_id = nan(1,n_stream);
            for i = 1 : n_stream
                sat(i) = rem(u_stream(i) ,1000);
                obs_id(i) = floor(u_stream(i)/1000);
                idx_res = this.obs_codes_id_obs == obs_id(i) & this.satellite_obs == sat(i);                
                if any(idx_res)
                    [~,idx_time] = ismember(this.time_obs.getEpoch(idx_res).getNominalTime.getRefTime(time_res.first.getMatlabTime),time_res.getNominalTime.getRefTime(time_res.first.getMatlabTime));
                    res_ph(idx_time, i) = this.res(idx_res);
                    res_id(idx_time, i) = find(idx_res);
                end
            end
            
        end
        
        function [res_pr, sat, obs_id] = getPrRes(this, rec)
            % get phase residual
            %
            % SYNTAX:  res_ph = getPhRes(this)
            if nargin <2
                rec = 1;
            end
            idx_rec = find(this.receiver_obs == rec);
            u_stream = unique(uint32(this.satellite_obs(idx_rec  & ~this.phase_obs )) + 1000*uint32(this.obs_codes_id_obs(idx_rec  & ~this.phase_obs )));
            n_stream = length(u_stream);
            time_res = this.time_obs.getNominalTime.getEpoch(idx_rec).minimum;
            duration = this.time_obs.getNominalTime.getEpoch(idx_rec).maximum - time_res;
            time_res.addSeconds(0:this.time_obs.getRate:duration);
            res_pr = nan(time_res.length,n_stream);
            sat = nan(1,n_stream);
            obs_id = nan(1,n_stream);
            for i = 1 : n_stream
                sat(i) = rem(u_stream(i) ,1000);
                obs_id(i) = floor(u_stream(i)/1000);
                idx_res = this.obs_codes_id_obs == obs_id(i) & this.satellite_obs == sat(i);
                if any(idx_res)
                    [~,idx_time] = ismember(this.time_obs.getEpoch(idx_res).getNominalTime.getRefTime(time_res.first.getMatlabTime),time_res.getNominalTime.getRefTime(time_res.first.getMatlabTime));
                    res_pr(idx_time,i) = this.res(idx_res);
                end
            end
            
        end
        
        function setPhFlag(this,rec,flag)
            % set phase outlier
            %
            % SYNTAX:  setPhFlag(this,rec,flag)
            idx_rec = find(this.receiver_obs == rec);
            u_stream = unique(uint32(this.satellite_obs(idx_rec  & this.phase_obs )) + 1000*uint32(this.obs_codes_id_obs(idx_rec  & this.phase_obs )));
            n_stream = length(u_stream);
            time_res = this.time.getNominalTime.getEpoch(idx_rec).minimum;
            duration = this.time.getNominalTime.getEpoch(idx_rec).maximum - time_res;
            time_res.addSeconds(0:this.time.getRate:duration);
            for i = 1 : n_stream
                sat(i) = rem(u_stream(i) ,1000);
                obs_id(i) = floor(u_stream(i)/1000);
                idx_res = find(this.obs_codes_id_obs == obs_id(i) & this.satellite_obs == sat(i));
                if any(idx_res)
                    [~,idx_time] = ismember(this.time_obs.getEpoch(idx_res).getNominalTime.getRefTime(time_res.first.getMatlabTime),time_res.getNominalTime.getRefTime(time_res.first.getMatlabTime));
                    this.outlier_obs(idx_res) = flag(idx_tim,i);
                end
            end
        end
        
        function setPrFlag(this,rec,flag)
            % set phase outlier
            %
            % SYNTAX:  setPhFlag(this,rec,flag)
            idx_rec = find(this.receiver_obs == rec);
            u_stream = unique(uint32(this.satellite_obs(idx_rec  & ~this.phase_obs )) + 1000*uint32(this.obs_codes_id_obs(idx_rec  & ~this.phase_obs )));
            n_stream = length(u_stream);
            time_res = this.time.getNominalTime.getEpoch(idx_rec).minimum;
            duration = this.time.getNominalTime.getEpoch(idx_rec).maximum - time_res;
            time_res.addSeconds(0:this.time.getRate:duration);
            for i = 1 : n_stream
                sat(i) = rem(u_stream(i) ,1000);
                obs_id(i) = floor(u_stream(i)/1000);
                idx_res = find(this.obs_codes_id_obs == obs_id(i) & this.satellite_obs == sat(i));
                if any(idx_res)
                    [~,idx_time] = ismember(this.time_obs.getEpoch(idx_res).getNominalTime.getRefTime(time_res.first.getMatlabTime),time_res.getNominalTime.getRefTime(time_res.first.getMatlabTime));
                    this.outlier_obs(idx_res) = flag(idx_tim,i);
                end
            end
        end
        
        
        function setUpSA(this, rec_work, id_sync, flag, param_selction, parametrization)
            % set up single point adjustment
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if strcmpi(flag,'???')
                o_tmp = rec_work.getObsSet('L??');
                o_tmp.keepEpochs(id_sync);
                this.addObsEq(rec_work, o_tmp, param_selction);
                o_tmp = rec_work.getObsSet('C??');
                o_tmp.keepEpochs(id_sync);
                this.addObsEq(rec_work, o_tmp, param_selction);
            else
                o_tmp = rec_work.getObsSet(flag);
                o_tmp.keepEpochs(id_sync);
                this.addObsEq(rec_work, o_tmp, param_selction);
            end
            if nargin < 6
                parametrization = LS_Parametrization();
            end
            this.unique_time = rec_work.time;
            %             ls_param.tropo(1) = LS_Parametrization.SPLINE_CUB;
            %             ls_param.tropo_opt = struct('spline_rate',900);
            %             ls_param.tropo_e(1) = LS_Parametrization.SPLINE_CUB;
            %             ls_param.tropo_e_opt = struct('spline_rate',3600);
            %             ls_param.tropo_n(1) = LS_Parametrization.SPLINE_CUB;
            %             ls_param.tropo_n_opt = struct('spline_rate',3600);
            this.bondParamsGenerateIdx(parametrization);
            
        end
        
        function setUpPPP(this, rec_work, id_sync, param_selction, parametrization)
            % set up precise point positionign
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if nargin < 4 || isempty(param_selction)
                param_selction = [this.PAR_REC_X;
                    this.PAR_REC_Y;
                    this.PAR_REC_Z;
                    this.PAR_REC_EB;
                    this.PAR_AMB;
                    this.PAR_REC_CLK;
                    this.PAR_TROPO;
                    this.PAR_TROPO_N;
                    this.PAR_TROPO_E;
                    this.PAR_IONO];  %
            end
            if nargin < 5
                parametrization = LS_Parametrization();
            end
            
            this.setUpSA(rec_work, id_sync, '???', param_selction, parametrization);
        end
        
        function setUpIonoFreePPP(this,rec_work,id_sync)
            % set up precise point positionign
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            param_selction = [this.PAR_REC_X;
                this.PAR_REC_Y;
                this.PAR_REC_Z;
                this.PAR_REC_EB;
                this.PAR_AMB;
                this.PAR_REC_CLK;
                this.PAR_TROPO;
                this.PAR_TROPO_N;
                this.PAR_TROPO_E];
            this.addObsEq(rec_work, rec_work.getPrefIonoFree('L','G'), param_selction);
            this.addObsEq(rec_work, rec_work.getPrefIonoFree('C','G'), param_selction);
            
            ls_param = LS_Parametrization();
            this.bondParamsGenerateIdx(ls_param);
        end
        
        function setUpNET(this, sta_list, coo_rate, flag, param_selction, parametrization)
            % set up single point adjustment
            %
            % SYNTAX:
            %   this.setUpSA(rec_work,id_sync,obs_type)
            if nargin < 4
                param_selction = [this.PAR_REC_X ;
                    this.PAR_REC_Y;
                    this.PAR_REC_Z;
                    this.PAR_REC_EB;
                    this.PAR_IONO;
                    %  this.PAR_REC_EB_LIN;
                    this.PAR_AMB;
                    this.PAR_REC_CLK;
                    this.PAR_TROPO;
                    this.PAR_TROPO_N;
                    this.PAR_TROPO_E;
                    this.PAR_SAT_CLK;
                    this.PAR_SAT_EB ];
            end
            if nargin < 5
                parametrization = LS_Parametrization();
            end
            % get time common at least to two receiver
            [p_time, id_sync] = Receiver_Work_Space.getSyncTimeExpanded(sta_list, coo_rate);
            id_rem = sum(~isnan(id_sync),2) <= 1;
            p_time.remEpoch(id_rem);
            id_sync(id_rem,:) = [];
            this.unique_time = p_time;
            
            
            % get observation types common to at least two reciver
            [o_codes, id_sync_o] = Receiver_Work_Space.getCommonObsCode(sta_list);
            id_rem_o = find(sum(~isnan(id_sync_o),2) <= 1)';
            %id_sync_o(id_rem_o,:) = nan;
            
            
            % add equations
            for r = 1 : length(sta_list)
                if strcmpi(flag,'???')
                    o_tmp = sta_list(r).work.getObsSet('L??');
                    o_tmp.keepEpochs(noNaN(id_sync(:,r)));
                    for o = id_rem_o
                        idx_rm = o_tmp.go_id == str2num(o_codes(o,4:6)) & strLineMatch(o_tmp.obs_code(:,2:4),o_codes(o,1:3));
                        if sum(idx_rm) > 0
                            o_tmp.removeColumn(idx_rm);
                            this.log.addMessage(sprintf('Observation %s from satellite %s is seen only from receiver %s : removing from network adjustement',o_codes(o,1:3),  Core.getConstellationCollector.getAntennaId(str2num(o_codes(o,4:6))), sta_list(r).getMarkerName4Ch));
                        end
                    end
                    if ~o_tmp.isEmpty
                        this.addObsEq(sta_list(r).work, o_tmp, param_selction);
                    end
                    o_tmp = sta_list(r).work.getObsSet('C??');
                    
                    o_tmp.keepEpochs(noNaN(id_sync(:,r)));
                    for o = id_rem_o
                        idx_rm = o_tmp.go_id == str2num(o_codes(o,4:6)) & strLineMatch(o_tmp.obs_code(:,2:4),o_codes(o,1:3));
                        if sum(idx_rm) > 0
                            o_tmp.removeColumn(idx_rm);
                            this.log.addMessage(sprintf('Observation %s from satellite %s is seen only from receiver %s : removing from network adjustement',o_codes(o,1:3),  Core.getConstellationCollector.getAntennaId(str2num(o_codes(o,4:6))), sta_list(r).getMarkerName4Ch));
                            
                        end
                    end
                    if ~o_tmp.isEmpty
                        this.addObsEq(sta_list(r).work, o_tmp, param_selction);
                    end
                else
                    o_tmp = sta_list(r).work.getObsSet(flag);
                    o_tmp.keepEpochs(noNaN(id_sync(:,r)));
                    for o = id_rem_o
                        idx_rm = o_tmp.go_id == str2num(o_codes(o,4:6)) & strLineMatch(o_tmp.obs_code(:,2:4),o_codes(o,1:3));
                        if sum(idx_rm) > 0
                            o_tmp.removeColumn(idx_rm);
                            this.log.addMessage(sprintf('Observation %s from satellite %s is seen only from receiver %s : removing from network adjustement',o_codes(o,1:3),  Core.getConstellationCollector.getAntennaId(str2num(o_codes(o,4:6))), sta_list(r).getMarkerName4Ch));
                            
                        end
                    end
                    this.addObsEq(sta_list(r).work, o_tmp, param_selction);
                end
                
            end
            this.bondParamsGenerateIdx(parametrization);
            this.absValRegularization(this.PAR_IONO, 1e4);
            
        end
        
        function [time_st, time_end] = getTimePar(this, idx)
            % get the parameter time as GPS_Time
            %
            % SYNTAX:
            %  [time_st, time_end] = getTimePar(this,idx)
            if nargin < 2
                idx = 1:size(this.time_par,1);
            end
            time_st = this.time_min.getCopy();
            time_st.addSeconds( double(this.time_par(idx,1)));
            time_end = this.time_min.getCopy();
            time_end.addSeconds( double(this.time_par(idx,2)));
        end
        
        
        function s0 = getSigma0Ph(this)
            % Get sigma0 of phase (PPP solution)
            %
            % SYNTAX:
            %  s0 = this.getSigma0Ph()
            s0 = mean(abs(this.res(this.phase_obs > 0)));
        end
    end
end
